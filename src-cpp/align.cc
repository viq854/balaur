//#include <bitpck.h>
//#include <varint.h>
//#include <vbyte.h>
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
void phase2_monolith(reads_t& reads, const ref_t& ref, std::vector<voting_results>& voting_results,  voting_stats& stats);
void finalize(reads_t& reads, const ref_t& ref, std::vector<voting_results>& results, voting_stats& stats, sam_writer_t& sam_io);
void free_kmer_buffers(std::vector<voting_task*>& encrypt_kmer_buffers);

void balaur_main(const char* fastaName, ref_t& ref, reads_t& reads, precomp_contig_io_t& contig_io, sam_writer_t& sam_io) {
	omp_set_num_threads(1);

	// --- phase 1 ---
	double start_time = omp_get_wtime();
	phase1_minhash(ref, reads);
	reads.free_minhash_data();

	if(params->load_mhi) {
		//std::cout << "Releasing index memory\n";
		//ref.index.release();
		if(params->precomp_contig_file_name.size() != 0) {
			contig_io.store_precomp_contigs(reads);
		}
	} else {
		contig_io.load_precomp_contigs(reads);
	}
	filter_candidate_contigs(reads);

	// --- phase 2 ---
	//load_kmer2_hashes(fastaName, ref, params);
	std::vector<voting_results> results;
	voting_stats stats;
	if(params->monolith) {
		phase2_monolith(reads, ref, results, stats);
	} else {
		//load_repeat_info(fastaName, ref, params);
		std::vector<voting_task*> encrypt_kmer_buffers;
		phase2_encryption(reads, ref, encrypt_kmer_buffers);
		
		phase2_voting(encrypt_kmer_buffers, results, stats);
		free_kmer_buffers(encrypt_kmer_buffers);
	}
	finalize(reads, ref, results, stats, sam_io);
	//eval(reads, ref);
	printf("****TOTAL ALIGNMENT TIME****: %.2f sec\n", omp_get_wtime() - start_time);	
}

// generate the read minhash firgerprints
// assemble candidate contigs
void phase1_minhash(const ref_t& ref, reads_t& reads) {
	printf("////////////// Phase 1: MinHash //////////////\n");
	double t = omp_get_wtime();
	///// ---- fingerprints ----
	//#pragma omp parallel for
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

	// compression (experimental)
	/*total_size = 0;
	for(size_t i = 0; i < encrypt_kmer_buffers.size(); i++) {
                voting_task* task = encrypt_kmer_buffers[i];
                kmer_cipher_t* data = task->get_data();
		int data_len = task->get_data_len();
		
		//uint32_t size = vbyte_compress_unsorted64(reinterpret_cast<uint64_t*>(data), out, data_len);
		std::vector<uint64_t> data_vec(data_len);
		std::vector<uint64_t> data_dec(data_len);
		for(int x = 0; x < data_len; x++) {
			data_vec[x] = data[x];
		}
		
		//size_t size = oroch::varint_codec<uint64_t>::space(data_vec.begin(), data_vec.end());
		size_t size = oroch::bitpck_codec<uint64_t>::space(data_len, 64);
		total_size += size;
		
		//std::vector<uint8_t> data_enc(size);
		//uint8_t *out = data_enc.data();
		//oroch::varint_codec<uint64_t>::encode(out, data_vec.begin(), data_vec.end());
		//const uint8_t* out2 = data_enc.data();
		//oroch::varint_codec<uint64_t>::decode(data_dec.begin(), data_dec.end(), out2);
	}
	printf("Total size: %.2f MB\n", ((float) total_size)/1024/1024);
	*/
}

void phase2_voting(std::vector<voting_task*>& encrypt_kmer_buffers, std::vector<voting_results>& results, voting_stats& stats) {
	printf("////////////// Phase 2: Voting //////////////\n");
	double t = omp_get_wtime();
	results.resize(encrypt_kmer_buffers.size());
	run_voting(encrypt_kmer_buffers, results, stats);
	printf("Total voting time: %.2f sec\n", omp_get_wtime() - t);

	// ---- determine the total communication size ----
	uint64 total_size = 0;
	for(size_t i = 0; i < results.size(); i++) {
		total_size += sizeof(results[i]); //TODO: more fine-grained
	} 
	printf("Total size: %.2f MB\n", ((float) total_size)/1024/1024);
}

void finalize(reads_t& reads, const ref_t& ref, std::vector<voting_results>& results, 
	voting_stats& stats, sam_writer_t& sam_writer) {
	printf("////////////// Finalize Mappings //////////////\n");
	double t = omp_get_wtime();
	for(size_t i = 0; i < results.size(); i++) {
		voting_results& task_out = results[i];
		read_t& r = reads.reads[task_out.rid];
		seq_t global_offset1 = 0;
		seq_t global_offset2 = 0;
		if(task_out.contig_id[0] >= 0) global_offset1 = r.ref_matches[task_out.contig_id[0]].pos; 
		if(task_out.contig_id[1] >= 0) global_offset2 = r.ref_matches[task_out.contig_id[1]].pos;
		task_out.convert2global_pos(global_offset1, global_offset2);

		if(task_out.kmer_matches.size() > 0) {	
			cipher_match_t m = task_out.kmer_matches[0];
			seq_t cpos = get_global_cipher_pos(m.ck.pos, m.ck.hash, global_offset1, ref.precomputed_kmer2_hashes, task_out.key1_xor_pad, task_out.key2_mult_pad, params->vanilla);
			seq_t rpos = m.rk.orig_pos;
			seq_t ref_pos = 0;
			if(cpos > rpos) {
				ref_pos = cpos - rpos;
			}
			task_out.global_pos[0] = ref_pos;
		}
		r.compare_and_update_best_aln(task_out.best_score, task_out.global_pos, task_out.rc);
	}

	int sum = 0;
	int n_nonzero = 0;
	//#pragma omp parallel for reduction(+:sum, n_nonzero)
	for(size_t i = 0; i < reads.reads.size(); i++) {
		read_t& r = reads.reads[i];
		if(r.top_aln.inlier_votes > 0) {
                        sum += r.top_aln.inlier_votes;
                        n_nonzero++;
                }
	}
	if(n_nonzero > 0) stats.avg_score = sum/n_nonzero;

	// ---- mapq score -----
	printf("Calculating the mapq scores...\n");
	//#pragma omp parallel for
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
			//enable for eval
			//if(r.top_aln.rc) {
			//	r.top_aln.ref_start += r.len;
			//}	
		}
	}
	
	// output the alignment results to SAM
	//store_alns_sam(reads, ref, params);
	sam_writer.write_sam_batch(reads, ref);
	printf("Total post-processing time: %.2f sec\n", omp_get_wtime() - t);
}

// allocate data transfer buffers
//storage for the ecrypted kmers
void allocate_encrypt_kmer_buffers(reads_t& reads, std::vector<voting_task*>& encrypt_kmer_buffers) {
	encrypt_kmer_buffers.reserve(reads.reads.size());
	int task_id = 0;
	for(size_t i = 0; i < reads.reads.size(); i++) {
		read_t& r = reads.reads[i];
		if(r.len <= 0) {
			std::cout << "ERROR: r.len \n";
			exit(-1);
		}

		if(!r.is_valid()) continue;
		if(r.n_match_f > 0) {
			const int n_batches = ceil(((float) r.n_match_f) / params->batch_size);
			for(int j = 0; j < n_batches; j++) {
					const int start = j * params->batch_size;
					const int end = (j == n_batches - 1) ? r.n_match_f : start + params->batch_size;
					voting_task* new_task = voting_task::alloc_voting_task(r.len, i, voting_task::strand_t::FWD, r.ref_matches, start, end);
					if(new_task != 0) {
						encrypt_kmer_buffers.push_back(new_task);
					}
			} 
		}
		if(r.ref_matches.size() - r.n_match_f > 0) {
			const int n_batches = ceil(((float) (r.ref_matches.size() - r.n_match_f)) / params->batch_size);
			for(int j = 0; j < n_batches; j++) {
					const int start =  r.n_match_f + j * params->batch_size;
					const int end = (j == n_batches - 1) ? r.ref_matches.size() : start + params->batch_size;
					voting_task* new_task = voting_task::alloc_voting_task(r.len, i, voting_task::strand_t::RC, r.ref_matches, start, end);
					if(new_task != 0) {
						encrypt_kmer_buffers.push_back(new_task);
					}
			} 
		}
	}
}

void free_kmer_buffers(std::vector<voting_task*>& encrypt_kmer_buffers) {
	//#pragma omp parallel for
	for(size_t i = 0; i < encrypt_kmer_buffers.size(); i++) {
		encrypt_kmer_buffers[i]->free();
	}
}

void populate_encrypt_kmer_buffers(reads_t& reads, const ref_t& ref, std::vector<voting_task*>& encrypt_kmer_buffers) {
	//#pragma omp parallel for (fix: rhashes access needs to be synchronized)
	for(size_t i = 0; i < encrypt_kmer_buffers.size(); i++) {
		voting_task* task = encrypt_kmer_buffers[i];
		read_t* r = &reads.reads[task->rid];
		const char* rseq;
		kmer_cipher_t** rhashes; 
		if(task->strand == voting_task::strand_t::FWD) {
			rseq =  r->seq.c_str();
			rhashes = &r->hashes_f;
		} else {
			rseq =  r->rc.c_str();
			rhashes = &r->hashes_rc;
		}
		if(params->vanilla) {
			generate_vanilla_ciphers(task->get_read(), rseq, r->len);
		} else {
			r->set_repeat_mask(params->k2, params->mask_repeat_nbrs ? params->k2 : 0);
			if(*rhashes == NULL) {
				generate_sha1_ciphers(task->get_read(), rseq, r->len, r->repeat_mask, task->strand);
				if(task->get_read_data_len() <= 0) {
					std::cout << "ERROR alloc read len <= 0 " << i << "\n";
					exit(-1);
				}
				*rhashes = new kmer_cipher_t[task->get_read_data_len()];
				memcpy(*rhashes, task->get_read(), sizeof(kmer_cipher_t)*task->get_read_data_len());
			} else {
				memcpy(task->get_read(), *rhashes, sizeof(kmer_cipher_t)*task->get_read_data_len());
			}
		}
		
		int contig_id = 0;
		for(int j = task->start; j < task->end; j++) {
			if(!r->ref_matches[j].valid) continue;
			if(params->vanilla) {
				 lookup_vanilla_ciphers(task->get_contig(contig_id), r->ref_matches[j].pos, r->ref_matches[j].len, ref.precomputed_kmer2_hashes);
			} else {
				lookup_sha1_ciphers(task->get_contig(contig_id), true, r->ref_matches[j].pos, r->ref_matches[j].len, ref.precomputed_kmer2_hashes, ref.precomputed_neighbor_repeats);
			}
			contig_id++;
#if(SIM_EVAL)
			r->get_sim_read_info(ref);
	 		if(pos_in_intv(r->ref_pos_r, r->ref_matches[j].pos, r->ref_matches[j].len) || pos_in_intv(r->ref_pos_l, r->ref_matches[j].pos, r->ref_matches[j].len))  {
				r->collected_true_hit = true;
				r->processed_true_hit = true;
				r->true_n_bucket_hits = r->ref_matches[j].n_diff_bucket_hits;
				task->true_cid = contig_id-1;
			}
#endif
		}
		// apply the task-specific keys
		if(!params->vanilla){
			task->key1_xor_pad = genrand64_int64();
        		task->key2_mult_pad = genrand64_int64();
			apply_keys(task->get_data(), task->get_data_len(), task->key1_xor_pad, task->key2_mult_pad);
		}
	}
}

// on-the-fly voting without buffering
void phase2_monolith(reads_t& reads, const ref_t& ref, std::vector<voting_results>& results, voting_stats& stats) {
	printf("////////////// Phase 2: MONOLITH //////////////\n");
	double t = omp_get_wtime();
	//#pragma omp parallel for
	for(size_t i = 0; i < reads.reads.size(); i++) {
		read_t& r = reads.reads[i];
		if(!r.is_valid()) continue;
		if(r.n_match_f > 0) {
			voting_task* task = voting_task::alloc_voting_task(r.len, i, voting_task::strand_t::FWD, r.ref_matches, 0, r.n_match_f);
			if(task != NULL) {
				generate_vanilla_ciphers(task->get_read(), r.seq.c_str(), r.len);
				int contig_id = 0;
				for(int j = 0; j < r.n_match_f; j++) {
					if(!r.ref_matches[j].valid) continue;
					lookup_vanilla_ciphers(task->get_contig(contig_id), r.ref_matches[j].pos, r.ref_matches[j].len, ref.precomputed_kmer2_hashes);
					contig_id++;
					
				#if(SIM_EVAL)
        		                r.get_sim_read_info(ref);
                        		if(pos_in_intv(r.ref_pos_r, r.ref_matches[j].pos, r.ref_matches[j].len) || pos_in_intv(r.ref_pos_l, r.ref_matches[j].pos, r.ref_matches[j].len))  {
                                		r.collected_true_hit = true;
                               			r.processed_true_hit = true;
                                		r.true_n_bucket_hits = r.ref_matches[j].n_diff_bucket_hits;
                                		task->true_cid = contig_id-1;
                        		}
				#endif		
	
				}
				voting_results res;
				res.rid = task->rid;
                		res.rc = task->strand;
				task->process(res);
				results.push_back(res);
				task->free();
			}
		}
		if(r.ref_matches.size() - r.n_match_f > 0) {
			voting_task* task = voting_task::alloc_voting_task(r.len, i, voting_task::strand_t::RC, r.ref_matches, r.n_match_f, r.ref_matches.size());
			if(task != NULL) {
				generate_vanilla_ciphers(task->get_read(), r.rc.c_str(), r.len);
				int contig_id = 0;
				for(int j = r.n_match_f; j < r.ref_matches.size(); j++) {
					if(!r.ref_matches[j].valid) continue;
					lookup_vanilla_ciphers(task->get_contig(contig_id), r.ref_matches[j].pos, r.ref_matches[j].len, ref.precomputed_kmer2_hashes);
					contig_id++;
				#if(SIM_EVAL)
                                        r.get_sim_read_info(ref);
                                        if(pos_in_intv(r.ref_pos_r, r.ref_matches[j].pos, r.ref_matches[j].len) || pos_in_intv(r.ref_pos_l, r.ref_matches[j].pos, r.ref_matches[j].len))  {
                                                r.collected_true_hit = true;
                                                r.processed_true_hit = true;
                                                r.true_n_bucket_hits = r.ref_matches[j].n_diff_bucket_hits;
                                                task->true_cid = contig_id-1;
                                        }
                                #endif				

				}
				voting_results res;
				res.rid = task->rid;
                                res.rc = task->strand;
                                task->process(res);
                                results.push_back(res);
				task->free();

				//if(res.n_true_votes > 0) {
				//	r.comp_votes_hit = res.n_true_votes;
				//}
			}
		}
	}
	printf("Total time: %.2f sec\n", omp_get_wtime() - t);
}

