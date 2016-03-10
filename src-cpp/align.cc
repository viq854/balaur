#include <float.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <utility>
#include <limits.h>
#include <queue>
#include <bitset>

#include "crypt.h"
#include "align.h"
#include "voting.h"
#include "contigs.h"
#include "index.h"
#include "io.h"
#include "hash.h"
#include "sam.h"
#include "lsh.h"

//////////// PRIVACY-PRESERVING READ ALIGNMENT ////////////
void phase1_minhash(const ref_t& ref, reads_t& reads);
void phase2_encryption(reads_t& reads, const ref_t& ref, std::vector<voting_task*>& encrypt_kmer_buffers);
void phase2_voting(std::vector<voting_task*>& encrypt_kmer_buffers, std::vector<voting_results>& results, voting_stats& stats);
void phase2_monolith(reads_t& reads, const ref_t& ref, std::vector<voting_results> voting_results,  voting_stats& stats);
void finalize(reads_t& reads, const ref_t& ref, std::vector<voting_results>& results, const voting_stats& stats);

void balaur_main(const char* fastaName, ref_t& ref, reads_t& reads) {
	// --- phase 1 ---
	double start_time = omp_get_wtime();
	phase1_minhash(ref, reads);
	
	if(params->load_mhi) {
		ref.index.release();
		if(params->precomp_contig_file_name.size() != 0) {
			store_precomp_contigs(params->precomp_contig_file_name.c_str(), reads);
		}
	} else {
		load_precomp_contigs(params->precomp_contig_file_name.c_str(), reads);
	}
	filter_candidate_contigs(reads);

	// --- phase 2 ---
	load_kmer2_hashes(fastaName, ref, params);
	load_repeat_info(fastaName, ref, params);
	//load_repeat_local(fastaName, ref, params);
	std::vector<voting_results> results;
	voting_stats stats;
	if(!params->monolith) {
		std::vector<voting_task*> encrypt_kmer_buffers;
		phase2_encryption(reads, ref, encrypt_kmer_buffers);
		phase2_voting(encrypt_kmer_buffers, results, stats);
	} else {
		phase2_monolith(reads, ref, results, stats);
	}
	finalize(reads, ref, results, stats);
	eval(reads, ref);
	printf("****TOTAL ALIGNMENT TIME****: %.2f sec\n", omp_get_wtime() - start_time);	
}

// generate the read minhash firgerprints
// assemble candidate contigs
void phase1_minhash(const ref_t& ref, reads_t& reads) {
	printf("////////////// Phase 1: MinHash //////////////\n");
	double t = omp_get_wtime();
	///// ---- fingerprints ----
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		r->minhashes_f.resize(params->h);
		r->minhashes_rc.resize(params->h);
		r->valid_minhash_f = minhash(r->seq, ref.high_freq_kmer_bitmap, r->minhashes_f);
		r->valid_minhash_rc = minhash(r->rc, ref.high_freq_kmer_bitmap, r->minhashes_rc);
	}
	printf("Runtime (fingerprints): %.2f sec\n", omp_get_wtime() - t);
	if(!params->load_mhi) return;
	
	// ---- candidate contigs ----
	assemble_candidate_contigs(ref, reads);
	printf("Runtime time (total): %.2f sec\n", omp_get_wtime() - t);
}

// encrypt the read and contig kmers
void allocate_encrypt_kmer_buffers(reads_t& reads, std::vector<voting_task*>& encrypt_kmer_buffers);
void populate_encrypt_kmer_buffers(reads_t& reads, const ref_t& ref, std::vector<voting_task*>& encrypt_kmer_buffers);
void phase2_encryption(reads_t& reads, const ref_t& ref, std::vector<voting_task*>& encrypt_kmer_buffers) {
	printf("////////////// Phase 2: Contig Encryption //////////////\n");
	double t1 = omp_get_wtime();
	allocate_encrypt_kmer_buffers(reads, encrypt_kmer_buffers);
	printf("Data alloc time: %.2f sec\n", omp_get_wtime() - t1);
	double t2 = omp_get_wtime();
	populate_encrypt_kmer_buffers(reads, ref, encrypt_kmer_buffers);
	printf("Encryption time: %.2f sec\n", omp_get_wtime() - t2);
	printf("Total time: %.2f sec\n", omp_get_wtime() - t1);
	
	// ---- determine the total communication size ----
	uint64 total_size = 0;
	uint64 total_contigs = 0;
	for(size_t i = 0; i < encrypt_kmer_buffers.size(); i++) {
		voting_task* task = encrypt_kmer_buffers[i];
		total_size += task->offsets[task->offsets.size()-1]*sizeof(kmer_cipher_t);
		total_contigs += task->offsets.size()-1;
	}
	printf("Total number of tasks: %lu \n", encrypt_kmer_buffers.size());
	printf("Total contigs: %llu \n", total_contigs);
	printf("Total size: %.2f MB\n", ((float) total_size)/1024/1024);
}

void phase2_voting(std::vector<voting_task*>& encrypt_kmer_buffers, std::vector<voting_results>& results, voting_stats& stats) {
	printf("////////////// Phase 2: Voting //////////////\n");
	double t = omp_get_wtime();
	run_voting(encrypt_kmer_buffers, results, stats);
	printf("Total voting time: %.2f sec\n", omp_get_wtime() - t);

	// ---- determine the total communication size ----
	uint64 total_size = 0;
	for(size_t i = 0; i < results.size(); i++) {
		total_size += sizeof(results[i]); //TODO: more fine-grained
	} 
	printf("Total size: %.2f MB\n", ((float) total_size)/1024/1024);
}

// on-the-fly voting without buffering
void phase2_monolith(reads_t& reads, const ref_t& ref, std::vector<voting_results> results,  voting_stats& stats) {
	printf("////////////// Phase 2: MONOLITH //////////////\n");
	double t = omp_get_wtime();
	int sum = 0;
	int n_nonzero = 0;
	for(size_t i = 0; i < reads.reads.size(); i++) {
		read_t& r = reads.reads[i];
		if(!r.is_valid()) continue;
		uint64 key1_xor_pad = genrand64_int64();
		uint64 key2_mult_pad = genrand64_int64();
		if(r.n_match_f > 0) {
			voting_task* new_task = voting_task::alloc_voting_task(r.len, i, voting_task::strand_t::FWD, r.ref_matches, 0, r.n_match_f);
			if(new_task != NULL) {
				generate_voting_kmer_ciphers_read(&new_task->data[1], r.seq.c_str(), r.len, key1_xor_pad, key2_mult_pad);
				for(int j = 0; j < r.n_match_f; j++) {
					if(!r.ref_matches[j].valid) continue;
					generate_voting_kmer_ciphers_ref(&new_task->data[1], ref.seq.c_str(), r.ref_matches[j].pos, r.ref_matches[j].len, key1_xor_pad, key2_mult_pad, ref);
				}
				voting_results res;
				new_task->process(res);
				results.push_back(res);
				if(res.best_score[0] > 0) {
					sum += res.best_score[0];
					n_nonzero++;
				}
			}
		}
		if(r.ref_matches.size() - r.n_match_f > 0) {
			voting_task* new_task = voting_task::alloc_voting_task(r.len, i, voting_task::strand_t::RC, r.ref_matches, r.n_match_f, r.ref_matches.size());
			if(new_task != NULL) {
				generate_voting_kmer_ciphers_read(&new_task->data[1], r.rc.c_str(), r.len, key1_xor_pad, key2_mult_pad);
				for(int j = r.n_match_f; j < r.ref_matches.size(); j++) {
					if(!r.ref_matches[j].valid) continue;
					generate_voting_kmer_ciphers_ref(&new_task->data[1], ref.seq.c_str(), r.ref_matches[j].pos, r.ref_matches[j].len, key1_xor_pad, key2_mult_pad, ref);
				}
			}
		}
	}
	printf("Total time: %.2f sec\n", omp_get_wtime() - t);
}

void finalize(reads_t& reads, const ref_t& ref, std::vector<voting_results>& results, const voting_stats& stats) {
	printf("////////////// Finalize Mappings //////////////\n");
	double t = omp_get_wtime();
	for(size_t i = 0; i < results.size(); i++) {
		voting_results& task_out = results[i];
		read_t& r = reads.reads[task_out.rid];
		task_out.convert2global_pos(r.ref_matches[task_out.contig_id[0]].pos, r.ref_matches[task_out.contig_id[1]].pos);
		r.compare_and_update_best_aln(task_out.best_score, task_out.global_pos, task_out.rc);
	}
	
	// ---- mapq score -----
	for(size_t i = 0; i < reads.reads.size(); i++) {
		read_t& r = reads.reads[i];
		r.top_aln.score = 0;
		// top > 0 and top != second best
		if(r.top_aln.inlier_votes > r.second_best_aln.inlier_votes) {
			// if sufficient votes were accumulated (lower thresholds for unique hit)
			if(r.top_aln.inlier_votes > params->votes_cutoff) {
				if(r.second_best_aln.inlier_votes < 0) r.second_best_aln.inlier_votes = 0;
				r.top_aln.score = params->mapq_scale_x*(r.top_aln.inlier_votes - r.second_best_aln.inlier_votes)/r.top_aln.inlier_votes;
				// scale by the distance from theoretical best possible votes
				if(stats.avg_score > 0 && params->enable_scale) {
					r.top_aln.score *= (float) r.top_aln.inlier_votes/(float) stats.avg_score;
				}
			}
			if(r.top_aln.rc) {
				r.top_aln.ref_start += r.len;
			}
		}
	}
	
	// output the alignment results to SAM
	store_alns_sam(reads, ref, params);
	printf("Total post-processing time: %.2f sec\n", omp_get_wtime() - t);
}

// allocate data transfer buffers
//storage for the ecrypted kmers
void allocate_encrypt_kmer_buffers(reads_t& reads, std::vector<voting_task*>& encrypt_kmer_buffers) {
	encrypt_kmer_buffers.reserve(reads.reads.size());
	for(size_t i = 0; i < reads.reads.size(); i++) {
		read_t& r = reads.reads[i];
		if(!r.is_valid()) continue;
		if(r.n_match_f > 0) {
			voting_task* new_task = voting_task::alloc_voting_task(r.len, i, voting_task::strand_t::FWD, r.ref_matches, 0, r.n_match_f);
			if(new_task != 0) {
				encrypt_kmer_buffers.push_back(new_task);
			}
		}
		if(r.ref_matches.size() - r.n_match_f > 0) {
			voting_task* new_task = voting_task::alloc_voting_task(r.len, i, voting_task::strand_t::RC, r.ref_matches, r.n_match_f, r.ref_matches.size());
			if(new_task != 0) {
				encrypt_kmer_buffers.push_back(new_task);
			}
		}
	}

}

void populate_encrypt_kmer_buffers(reads_t& reads, const ref_t& ref, std::vector<voting_task*>& encrypt_kmer_buffers) {
	for(size_t i = 0; i < encrypt_kmer_buffers.size(); i++) {
		voting_task* task = encrypt_kmer_buffers[i];
		read_t& r = reads.reads[task->rid];
		uint64 key1_xor_pad = genrand64_int64();
		uint64 key2_mult_pad = genrand64_int64();
		if(task->strand == voting_task::strand_t::FWD) {
			generate_voting_kmer_ciphers_read(task->get_read(), r.seq.c_str(), r.len, key1_xor_pad, key2_mult_pad);
		} else {
			generate_voting_kmer_ciphers_read(task->get_read(), r.rc.c_str(), r.len, key1_xor_pad, key2_mult_pad);
		}
		for(int j = task->start; j < task->end; j++) {
			if(!r.ref_matches[j].valid) continue;
			int contig_id = j - task->start;
			generate_voting_kmer_ciphers_ref(task->get_contig(j - task->start), ref.seq.c_str(), r.ref_matches[j].pos, r.ref_matches[j].len, key1_xor_pad, key2_mult_pad, ref);
		}
	}
}
