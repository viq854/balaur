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
#include <assert.h>
#include <queue>
#include <unordered_map>

#include "index.h"
#include "io.h"
#include "hash.h"
#include "cluster.h"
#include <bitset>

#include "ksw.h"

int eval_read_hit(ref_t& ref, read_t* r, const index_params_t* params);

struct heap_entry_t {
	seq_t pos;
	uint32 tid;
	uint32 next_idx;
};

void heap_sort(heap_entry_t* heap, int n) {
	heap_entry_t tmp;
	int i, j;
	for(j = 1; j < n; j++) {
		tmp = heap[j];
		for(i = j - 1; (i >= 0) && (heap[i].pos > tmp.pos); i--) {
			heap[i+1] = heap[i];
		}
		heap[i+1] = tmp;
	}
}

void heap_update(heap_entry_t* heap, uint32 n) {
	uint32 i = 0;
	uint32 k = i;
	heap_entry_t tmp = heap[i];
	while((k = (k << 1) + 1) < n) {
		if(k != (n - 1) && (heap[k].pos >= heap[k+1].pos)) ++k;
		if(heap[k].pos >= tmp.pos) break;
		heap[i] = heap[k];
		i = k;
	}
	heap[i] = tmp;
}

// min-heap

void swap(heap_entry_t*x, heap_entry_t*y) {
	heap_entry_t temp = *x;  *x = *y;  *y = temp;
}

void sift_down(heap_entry_t* heap, uint32 n, int i) {
	int l = 2*i + 1;
	int r = 2*i + 2;
	int min = i;
	if (l < n && heap[l].pos < heap[i].pos) {
		min = l;
	}
	if (r < n && heap[r].pos < heap[min].pos) {
		min = r;
	}
	if (min != i) {
		swap(&heap[i], &heap[min]);
	    sift_down(smallest);
	}
}

void heap_create(heap_entry_t* heap, int n) {
	int i = (n-1)/2;
	while(i >= 0) {
		sift_down(heap, n, i);
		i--;
	}
}


inline void heap_update_memmove(heap_entry_t* heap, uint32 n) {
	if(n <= 1) return;
	if(heap[1].pos >= heap[0].pos) return; // nothing to shift

	heap_entry_t tmp = heap[0];
	if(heap[n-1].pos <= heap[0].pos) { // shift all
		memmove(heap, heap+1, (n-1)*sizeof(heap_entry_t));
		heap[n-1] = heap[0];
	} else { // binary search: find the first element > tmp
		int i, j, k;
		i = 1;
		j = n;
		while(i != j) {
			k = (i + j)/2;
			if (heap[k].pos < tmp.pos) {
				i = k + 1;
			} else {
				j = k;
			}
		}
		memmove(heap, heap+1, (i)*sizeof(heap_entry_t));
		heap[i-1] = tmp;
	}
}

// output matches (ordered by the number of projections matched)
void collect_read_hits_contigs_inssort_pqueue(ref_t& ref, read_t* r, const index_params_t* params) {
	r->ref_matches.resize(params->n_tables);
	int n_best_hits = 0; // best number of table hits found so far

	// priority heap of matched positions
	heap_entry_t heap[params->n_tables];
	uint32 heap_size = 0;

	// push the first entries in each sorted bucket onto the heap
	for(uint32 t = 0; t < params->n_tables; t++) { // for each table
		if(r->ref_bucket_matches_by_table[t] == NULL) {
			continue;
		} else {
			heap[heap_size].pos = (*r->ref_bucket_matches_by_table[t])[0];
			heap[heap_size].tid = t;
			heap[heap_size].next_idx = 1;
		}
		heap_size++;
	}
	if(heap_size == 0) return; // all the matched buckets are empty
	//heap_sort(heap, heap_size); // build heap
	heap_create(heap, heap_size);

	int n_diff_table_hits = 0;
	uint32 len = 0;
	uint last_pos = -1;
	std::bitset<32> occ;
	while(heap_size > 0) {
		heap_entry_t e = heap[0]; // get min
		if(last_pos == (uint32) -1 || (e.pos <= last_pos + params->contig_gap)) { // first contig or extending contig
			if(!occ.test(e.tid)) {
				n_diff_table_hits++;
			}
			len += e.pos - last_pos;
			last_pos = e.pos;
			occ.set(e.tid);
		} else { // found a boundary, store last contig
			if(n_diff_table_hits >= (int) params->min_n_hits && n_diff_table_hits >= (n_best_hits - 1)) {
				if(n_diff_table_hits > n_best_hits) { // if more hits than best so far
					n_best_hits = n_diff_table_hits;
				}
				if(r->ref_matches[n_diff_table_hits-1].size() < params->max_best_hits) {
					ref_match_t rm(last_pos, len);
					r->ref_matches[n_diff_table_hits-1].push_back(rm);
				}
			}
			// start a new contig
			n_diff_table_hits = 1;
			len = 0;
			last_pos = e.pos;
			occ.reset();
			occ.set(e.tid);
		}
		// push the next match from this bucket
		if(e.next_idx < (*r->ref_bucket_matches_by_table[e.tid]).size()) {
			heap[0].pos = (*r->ref_bucket_matches_by_table[e.tid])[e.next_idx];
			heap[0].tid = e.tid;
			heap[0].next_idx = e.next_idx+1;
			//heap_update(heap, heap_size);
			sift_down(heap, heap_size, 0);
		} else { // no more entries in this bucket
			heap[0].pos = UINT_MAX;
			//heap_update(heap, heap_size);
			sift_down(heap, heap_size, 0);
			heap_size--;
		}
	}

	// add the last position
	if(last_pos != (uint32) -1 && n_diff_table_hits >= (int) params->min_n_hits && n_diff_table_hits >= (n_best_hits - 1)) {
		if(n_diff_table_hits > n_best_hits) { // if more hits than best so far
			n_best_hits = n_diff_table_hits;
		}
		if(r->ref_matches[n_diff_table_hits-1].size() < params->max_best_hits) {
			ref_match_t rm(last_pos, len);
			r->ref_matches[n_diff_table_hits-1].push_back(rm);
		}
	}

	r->best_n_hits = (n_best_hits > 0) ? n_best_hits - 1 : 0;
	assert(r->best_n_hits < params->n_tables);
	std::vector< VectorSeqPos* >().swap(r->ref_bucket_matches_by_table); //release memory
}

struct seed_t {
	seq_t ref_pos;
	seq_t read_pos;
	uint32 len;
	uint32 ref_contig_idx;
	seed_t(seq_t _ref_pos, seq_t _read_pos, uint32 _len, uint32 _ref_contig_idx) : ref_pos(_ref_pos), read_pos(_read_pos), len(_len), ref_contig_idx(_ref_contig_idx) {}
};

typedef std::vector<seed_t> SeedChain;

uint32 get_chain_weight(const SeedChain& seeds) {
	int32_t chain_end_pos = 0;
	uint32 weight_read = 0;
	for (uint32 i = 0; i < seeds.size(); i++) {
		const seed_t s = seeds[i];
		if (s.read_pos >= chain_end_pos) {
			weight_read += s.len;
		} else if (s.read_pos + s.len > chain_end_pos) {
			weight_read += s.read_pos + s.len - chain_end_pos;
		}
		chain_end_pos = chain_end_pos > s.read_pos + s.len ? chain_end_pos : s.read_pos + s.len;
	}
	uint32 weight_ref = 0;
	chain_end_pos = 0;
	for (uint32 i = 0; i < seeds.size(); i++) {
		const seed_t s = seeds[i];
		if (s.ref_pos >= chain_end_pos) {
			weight_read += s.len;
		} else if (s.ref_pos + s.len > chain_end_pos) {
			weight_read += s.ref_pos + s.len - chain_end_pos;
		}
		chain_end_pos = chain_end_pos > s.ref_pos + s.len ? chain_end_pos : s.ref_pos + s.len;
	}
	return weight_ref > weight_read ? weight_ref : weight_read;
}

// return 1 if the seed is merged into the chain
static bool add_seed(std::vector<seed_t>& seeds, const seed_t s, const index_params_t* params) {
	uint32 read_end, ref_end, x, y;
	const seed_t last = seeds[seeds.size()-1];
	read_end = last.read_pos + last.len;
	ref_end = last.ref_pos + last.len;
	if (s.read_pos >= seeds[0].read_pos && s.read_pos + s.len <= read_end && s.ref_pos >= seeds[0].ref_pos && s.ref_pos + s.len <= ref_end)
		return true; // seed already in chain

	x = s.read_pos - last.read_pos;
	y = s.ref_pos - last.ref_pos;
	if (x >= 0 && x - y <= params->bandw && y - x <= params->bandw && x - last.len < params->max_chain_gap && y - last.len < params->max_chain_gap) { // grow the chain
		seeds.push_back(s);
		return true;
	}
	return false;
}

struct comp_seeds
{
    bool operator()(const seed_t& a, const seed_t& b) const {
    	return a.len > b.len;
    }
};

struct comp_chains
{
    bool operator()(const SeedChain& a, const SeedChain& b) const {
    	return get_chain_weight(a) > get_chain_weight(b);
    }
};

#define CONTIG_PADDING 100

bool seed2alignment(const seed_t s, const ref_t& ref, read_t* r, const index_params_t* params) {
	aln_t* aln = &r->aln;
	ref_match_t ref_contig = r->ref_matches[r->best_n_hits][s.ref_contig_idx];
	int hit_len = r->len + ref_contig.len;
	seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;
	seq_t padded_hit_offset = (hit_offset >= CONTIG_PADDING) ? hit_offset - CONTIG_PADDING : 0;

	//printf("HIT: contig_len = %u contig_pos = %u offset = %u, len = %u \n", ref_contig.len, ref_contig.pos, hit_offset, hit_len);
	//printf("SEED: seed_len = %u seed_read_pos = %u seed_ref_pos = %u\n", s->len, s->read_pos, s->ref_pos);
	if (s.read_pos > 0) { // left extension
		int qle, tle, gtle, gscore;

		// populate the ref and read pieces
		uint32_t ref_len = s.ref_pos - padded_hit_offset;
		uint8_t ref_seq[ref_len];
		uint8_t read_seq[s.read_pos];
		for (uint32 i = 0; i < s.read_pos; i++) {
			read_seq[i] = r->seq[s.read_pos - 1 - i];
		}
		for (uint32 i = 0; i < ref_len; i++) {
			ref_seq[i] = ref.seq[padded_hit_offset + ref_len - 1 - i];
		}
		//printf("LEFT EXT: read_pos %u ref_pos %u ref_len %u\n", s->read_pos, s->ref_pos, ref_len);
		int offset;
		aln->score = ksw_extend2(s.read_pos, read_seq, ref_len, ref_seq, 5, params->score_matrix, params->gap_open, params->gap_extend, params->gap_open, params->gap_extend, params->bandw,
					0, params->zdrop, s.len*params->match, &qle, &tle, &gtle, &gscore, &offset);

		if (gscore >= 0) { // to-end
			aln->read_start = 0;
			aln->ref_start = s.ref_pos - gtle;
			aln->truesc = gscore;
		} else {
			// TODO: did not reach the end of the query!
			return false;
		}
	} else {
		aln->score = s.len * params->match;
		aln->read_start = 0;
		aln->ref_start = s.ref_pos;
	}

	if (s.read_pos + s.len != r->len) { // right extension
		int qle, tle, gtle, gscore;
		int sc0 = aln->score;
		uint32 read_len = r->len - (s.read_pos + s.len);
		uint32 ref_len = hit_offset + hit_len - (s.ref_pos + s.len);
		//printf("RIGHT EXT: read_pos %u read_len %u ref_pos %u ref_len %u\n", s.read_pos, read_len, s.ref_pos, ref_len);
		int offset;
		aln->score = ksw_extend2(read_len, (const unsigned char*) &(r->seq.c_str()[s.read_pos + s.len]), ref_len, (const unsigned char*)&(ref.seq.c_str()[s.ref_pos + s.len]), 5, params->score_matrix, params->gap_open, params->gap_extend, params->gap_open, params->gap_extend, params->bandw,
				0, params->zdrop, sc0, &qle, &tle, &gtle, &gscore, &offset);

		/*printf("REF: \n");
		for(uint32 i = 0; i < ref_len + s.len; i++) {
			printf("%c", iupacChar[(int)ref.seq.c_str()[s.ref_pos + i]]);
	   	} printf("\n");

	    printf("READ: \n");
		for(uint32 i = 0; i < read_len + s.len; i++) {
			printf("%c", iupacChar[(int)r->seq.c_str()[s.read_pos + i]]);
		} printf("\n");
		*/

		// similar to the above
		if(gscore >= 0) { // to-end extension
			aln->read_end = r->len;
			aln->ref_end = s.ref_pos + s.len + gtle;
			aln->truesc += gscore - sc0;
		} else {
			// TODO: did not reach the end of the query!
			return false;
		}
	} else {
		aln->read_end = r->len;
		aln->ref_end = s.ref_pos + s.len;
	}
	return true;
}

void global_alignment(const ref_match_t ref_contig, const ref_t& ref, read_t* r, const index_params_t* params) {
	aln_t* aln = &r->aln;
	int hit_len = r->len + ref_contig.len;
	seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;
	seq_t padded_hit_offset = (hit_offset >= CONTIG_PADDING) ? hit_offset - CONTIG_PADDING : 0;
	uint32_t ref_len = hit_len + 2*CONTIG_PADDING;

	int n_cigar;
	uint32_t* cigar;

	aln->score = ksw_global(r->len,  (const unsigned char*) r->seq.c_str(),
			ref_len, (const unsigned char*)ref.seq.c_str(),
			5, params->score_matrix, params->gap_open, params->gap_extend, params->bandw,
			&n_cigar, &cigar);
}

#define MAX_SEED_HITS 10

void process_read_hits_global(ref_t& ref, read_t* r, const index_params_t* params) {

	for(uint32 i = 0; i < r->ref_matches[r->best_n_hits].size(); i++) {
		ref_match_t ref_contig = r->ref_matches[r->best_n_hits][i];
		global_alignment(ref_contig, ref, r, params);
	}
}

void process_read_hits_se(ref_t& ref, read_t* r, const index_params_t* params) {
	// 1. index the read sequence: generate and store all kmers
	std::unordered_map<std::string, std::vector<uint32>> kmer2pos;
	for(uint32 i = 0; i < (r->len - params->k + 1); i++) {
		std::vector<uint32>& pos_vec = kmer2pos[std::string(&r->seq[i], params->k)];
		if(pos_vec.size() == 0) pos_vec.reserve(MAX_SEED_HITS);
		pos_vec.push_back(i);
	}

	// 2. find seeds shared between the reference and the read and assemble seed chains
	std::vector<SeedChain> chains;
	chains.reserve(r->ref_matches[r->best_n_hits].size());
	for(uint32 i = 0; i < r->ref_matches[r->best_n_hits].size(); i++) {
		// REF CANDIDATE CONTIG
		ref_match_t ref_contig = r->ref_matches[r->best_n_hits][i];
		int hit_len = r->len + ref_contig.len;
		seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;

		SeedChain* last_chain = NULL;
		uint32 step = 1;
		for(uint32 j = 0; j <= hit_len - params->k; j += step) { // REF KMER
			std::vector<uint32>& ref_kmer_hits = kmer2pos[std::string(&ref.seq[hit_offset + j], params->k)];
			if(ref_kmer_hits.size() == 0 || ref_kmer_hits.size() > MAX_SEED_HITS) continue;  //skip if absent or too frequent
			uint32 min_step = UINT_MAX;
			for(uint32 read_pos_idx = 0; read_pos_idx < ref_kmer_hits.size(); read_pos_idx++) { // READ HIT POS
				uint32 read_pos = ref_kmer_hits[read_pos_idx];
				uint32 e;
				for(e = params->k; e < (hit_len-j); e++) { // extend as much as possible to the right
					if(read_pos + e >= r->len) break; // reached the end of the read
					if(ref.seq[hit_offset + j + e] != r->seq[read_pos + e]) break; // reached a mismatch
				}

				seed_t s(hit_offset + j, read_pos, e, i); // new seed
				if(last_chain == 0 || !add_seed(*last_chain, s, params)) { // check if seed should be added to the last chain or if it should be a new chain
					chains.push_back(SeedChain());
					chains[chains.size()-1].reserve(10);
					chains[chains.size()-1].push_back(s);
					last_chain = &chains[chains.size()-1];
				}

				if(e < min_step) {
					min_step = e;
				}
			}
			// skip positions guaranteed to mismatch with the read
			step = min_step != UINT_MAX ? min_step : 1;
		}
	}
	if(chains.size() == 0) return;

	// 3. extend longest seeds with longest chains
	std::sort(chains.begin(), chains.end(), comp_chains());
	bool matched = false;
	for(uint32 i = 0; i < chains.size(); i++) {
		//printf("chain %u len %u \n", i, chains[i].size());
		std::sort(chains[i].begin(), chains[i].end(), comp_seeds());
		for(uint32 j = 0; j < chains[i].size(); j++) {
			//printf("seed %u len %u \n", j, chains[i][j].len);
			if(seed2alignment(chains[i][j], ref, r, params)) {
				matched = true;
				break;
			}
		}
		if(matched) break;
	}

	/*for (uint32 i = 0; i < chains.size(); i++) {
		chain_t p = chains[i];
		printf("* Found CHAIN(%d): num_seeds=%u \n", i, p.seeds.size());
		for (uint32 j = 0; j < p.seeds.size(); j++) {
			uint32 pos;
			printf("\t seed_len=%u;read_pos=%u;ref_pos=%u", p.seeds[j].len,  p.seeds[j].read_pos, p.seeds[j].ref_pos);
		}
		printf("\n");
	}*/
}

void process_read_hits_se_opt(ref_t& ref, read_t* r, const index_params_t* params) {
	// 1. index the read sequence: generate and store all kmers
	std::unordered_map<std::string, uint32> kmer2pos;
	kmer2pos[std::string(&r->seq[0], params->k)] = (uint32) -1;
	for(uint32 i = 1; i < (r->len - params->k + 1); i++) {
		uint32& pos = kmer2pos[std::string(&r->seq[i], params->k)];
		if(pos == 0) pos = i; // save the first occurrence
	}

	// 2. find seeds shared between the reference and the read
	bool matched = false;
	for(uint32 i = 0; i < r->ref_matches[r->best_n_hits].size(); i++) { // REF CANDIDATE CONTIG
		ref_match_t ref_contig = r->ref_matches[r->best_n_hits][i];
		seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;
		seq_t padded_hit_offset = (hit_offset >= CONTIG_PADDING) ? hit_offset - CONTIG_PADDING : 0;
		uint32 search_len = ref_contig.len + 2*CONTIG_PADDING + 200;

		uint32 step = 1;
		for(uint32 j = 0; j <= search_len; j += step) { // REF KMER
			step = 1;
			uint32& read_pos = kmer2pos[std::string(&ref.seq[padded_hit_offset + j], params->k)];
			if(read_pos == 0) {  //skip if absent
				continue;
			}
			if(read_pos == (uint32) -1) {
				read_pos = 0;
			}
			uint32 e;
			for(e = params->k; e < (search_len-j); e++) { // extend as much as possible to the right
				if(read_pos + e >= r->len) break; // reached the end of the read
				if(ref.seq[padded_hit_offset + j + e] != r->seq[read_pos + e]) break; // reached a mismatch
			}
			seed_t s(padded_hit_offset + j, read_pos, e, i); // new seed
			// 3. extend the seed
			if(seed2alignment(s, ref, r, params)) {
				matched = true;
				break;
			}
			step = e; // skip positions guaranteed to mismatch with the read
		}
		if(matched) break;
	}
}

#define MAX_TOP_HITS 100

void align_reads_minhash(ref_t& ref, reads_t& reads, const index_params_t* params) {
	printf("**** SRX Alignment: MinHash ****\n");

	uint32 max_windows_matched = 0;
	uint32 total_windows_matched = 0;
	uint32 total_top_contigs = 0;
	uint32 total_contigs_length = 0;
	uint32 diff_num_top_hits = 0;
	double start_time = omp_get_wtime();
	omp_set_num_threads(params->n_threads); // split the reads across the threads
	#pragma omp parallel reduction(+:total_windows_matched, total_top_contigs, diff_num_top_hits, total_contigs_length) reduction(max:max_windows_matched)
	{
		int tid = omp_get_thread_num();
		int n_threads = omp_get_num_threads();
		seq_t chunk_start = (reads.reads.size()/n_threads)*tid;
		seq_t chunk_end = (reads.reads.size()/n_threads)*(tid + 1);
		if(tid == n_threads - 1) {
			chunk_end = reads.reads.size();
		}
		printf("Thread %d range: %u %u \n", tid, chunk_start, chunk_end);

		for (uint32 i = chunk_start; i < chunk_end; i++) { // for each read of the thread's chunk
			if((i - chunk_start) % 10000 == 0 && (i - chunk_start) != 0) {
				printf("Thread %d processed %u reads \n", tid, i - chunk_start);
			}
			read_t* r = &reads.reads[i];
			r->best_n_hits = 0;
			if(!r->valid_minhash) continue;
			//r->ref_bucket_id_matches_by_table.resize(params->n_tables);
			r->ref_bucket_matches_by_table.resize(params->n_tables);
			for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
				minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes, params->sketch_proj_indices, t*params->sketch_proj_len);
				uint32 bucket_index = ref.hash_tables[t].bucket_indices[bucket_hash];
				if(bucket_index == ref.hash_tables[t].n_buckets) {
					continue; // no reference window fell into this bucket
				}
				r->ref_bucket_matches_by_table[t] = &ref.hash_tables[t].buckets_data_vectors[bucket_index];
			}
			collect_read_hits_contigs_inssort_pqueue(ref, r, params);

			if(r->best_n_hits > 0 && r->ref_matches[r->best_n_hits].size() < MAX_TOP_HITS) {
				//process_read_hits_se_opt(ref, r, params);
				//process_read_hits_global(ref, r, params);
			}

			// stats
			uint32 n_contigs = 0;
			uint32 top_contigs = 0;
			uint32 contig_length = 0;
			for(uint32 t = 0; t < params->n_tables; t++) {
				n_contigs += r->ref_matches[t].size();
				if(t == r->best_n_hits) {
					top_contigs += n_contigs;;
				}
				for(uint32 j = 0; j < r->ref_matches[t].size(); j++) {
				   ref_match_t match = r->ref_matches[t][j];
				   contig_length += match.len;
			   }
			}
			if(n_contigs > max_windows_matched) {
				max_windows_matched = n_contigs;
			}
			total_windows_matched += n_contigs;
			total_top_contigs += top_contigs;
			total_contigs_length += contig_length;
		}
	}
	double end_time = omp_get_wtime();
	printf("Collected read hits \n");

	printf("Evaluating read hits... \n");
	int valid_hash = 0;
	int acc_hits = 0;
	int acc_top = 0;
	int acc_dp = 0;
	//#pragma omp parallel for reduction(+:valid_hash, acc_hits, acc_top)
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		if(!reads.reads[i].valid_minhash) continue;
		eval_read_hit(ref, &reads.reads[i], params);
		acc_hits += reads.reads[i].acc;
		acc_top += reads.reads[i].top_hit_acc;
		acc_dp += reads.reads[i].dp_hit_acc;
		valid_hash += reads.reads[i].valid_minhash;
	}

	printf("Max number of windows matched by read %u \n", max_windows_matched);
	printf("Avg number of windows matched per read %.8f \n", (float) total_windows_matched/acc_hits);
	printf("Avg number of top contigs matched per read %.8f \n", (float) total_top_contigs/acc_hits);
	printf("Avg contig length per read %.8f \n", (float) total_contigs_length/total_windows_matched);
	printf("Avg diff of top 2 hits per read %.8f \n", (float) diff_num_top_hits/reads.reads.size());
	printf("Total number of accurate hits matching top = %d \n", acc_top);
	printf("Total number of accurate hits found = %d \n", acc_hits);
	printf("Total DP number of accurate hits found = %d \n", acc_dp);
	printf("Total search time: %f sec\n", end_time - start_time);

}

int eval_read_hit(ref_t& ref, read_t* r, const index_params_t* params) {
   unsigned int seq_id, pos_l, pos_r;
   int strand;
   parse_read_mapping(r->name.c_str(), &seq_id, &pos_l, &pos_r, &strand);
   seq_id = seq_id - 1;
   r->top_hit_acc = 0;

   assert(seq_id >= 0);
   if(ref.subsequence_offsets.size() > 1) {
	   assert(seq_id < ref.subsequence_offsets.size());
	   pos_l += ref.subsequence_offsets[seq_id]; // convert to global id
   }

   assert(r->best_n_hits < params->n_tables);
   for(int i = r->best_n_hits; i >= 0; i--) {
	   for(uint32 j = 0; j < r->ref_matches[i].size(); j++) {
		   ref_match_t match = r->ref_matches[i][j];
		   if(pos_l >= match.pos - match.len - 130 && pos_l <= match.pos + 130) {
			   r->acc = 1;
			   break;
		   }
	   }
	   if(r->acc == 1) {
		   if((uint32) i == r->best_n_hits) {
			  r->top_hit_acc = 1;
		   }
		   break;
	   }
   }

   if(pos_l >= r->aln.ref_start - 30 && pos_l <= r->aln.ref_start + 30) {
	   r->dp_hit_acc = 1;
   }

   return (r->acc == 1);
}

