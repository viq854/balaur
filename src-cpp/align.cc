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
		   if(pos_l >= match.pos - match.len - 30 && pos_l <= match.pos + 30) {
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
   return (r->acc == 1);
}

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
	heap_sort(heap, heap_size);

	int n_diff_table_hits = 0;
	uint32 len = 0;
	uint last_pos = -1;
	std::bitset<32> occ;
	while(heap_size > 0) {
		heap_entry_t e = heap[0];
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
			heap_update(heap, heap_size);
		} else { // no more entries in this bucket
			heap[0].pos = UINT_MAX;
			heap_update(heap, heap_size);
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

struct chain_t {
	int32_t pos;
	std::vector<seed_t> seeds;
};

uint32 get_chain_weight(const chain_t* c) {
	int32_t chain_end_pos = 0;
	uint32 weight_read = 0;
	for (uint32 i = 0; i < c->seeds.size(); i++) {
		const seed_t s = c->seeds[i];
		if (s.read_pos >= chain_end_pos) {
			weight_read += s.len;
		} else if (s.read_pos + s.len > chain_end_pos) {
			weight_read += s.read_pos + s.len - chain_end_pos;
		}
		chain_end_pos = chain_end_pos > s.read_pos + s.len ? chain_end_pos : s.read_pos + s.len;
	}
	uint32 weight_ref = 0;
	chain_end_pos = 0;
	for (uint32 i = 0; i < c->seeds.size(); i++) {
		const seed_t s = c->seeds[i];
		if (s.ref_pos >= chain_end_pos) {
			weight_read += s.len;
		} else if (s.ref_pos + s.len > chain_end_pos) {
			weight_read += s.ref_pos + s.len - chain_end_pos;
		}
		chain_end_pos = chain_end_pos > s.ref_pos + s.len ? chain_end_pos : s.ref_pos + s.len;
	}
	return weight_ref > weight_read ? weight_ref : weight_read;
}


#include "kbtree.h"

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, chain_t, chain_cmp)

// return 1 if the seed is merged into the chain
static int test_and_merge(const index_params_t* params, chain_t *c, const seed_t *p) {
	uint32 read_end, ref_end, x, y;
	const seed_t *last = &c->seeds[c->seeds.size()-1];
	read_end = last->read_pos + last->len;
	ref_end = last->ref_pos + last->len;
	if (p->read_pos >= c->seeds[0].read_pos && p->read_pos + p->len <= read_end && p->ref_pos >= c->seeds[0].ref_pos && p->ref_pos + p->len <= ref_end)
		return 1; // contained seed; do nothing

	x = p->read_pos - last->read_pos; // always non-negtive
	y = p->ref_pos - last->ref_pos;
	if (y >= 0 && x - y <= params->bandw && y - x <= params->bandw && x - last->len < params->max_chain_gap && y - last->len < params->max_chain_gap) { // grow the chain
		c->seeds.push_back(*p);
		return 1;
	}
	return 0; // request to add a new chain
}

void seed2alignment(seed_t* s, ref_t& ref, read_t* r, const index_params_t* params) {
	aln_t* aln = &r->aln;
	ref_match_t ref_contig = r->ref_matches[r->best_n_hits][s->ref_contig_idx];
	int hit_len = r->len + ref_contig.len;
	seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;

	if (s->read_pos > 0) { // left extension
		int qle, tle, gtle, gscore;

		// populate the ref and read pieces
		uint32_t ref_len = s->ref_pos - hit_offset;
		uint8_t ref_seq[ref_len];
		uint8_t read_seq[s->read_pos];
		for (uint32 i = 0; i < s->read_pos; i++) {
			read_seq[i] = r->seq[s->read_pos - 1 - i];
		}
		for (uint32 i = 0; i < ref_len; i++) {
			ref_seq[i] = ref.seq[hit_offset + ref_len - 1 - i];
		}
		int offset;
		aln->score = ksw_extend2(s->read_pos, read_seq, ref_len, ref_seq, 5, params->score_matrix, params->gap_open, params->gap_extend, params->gap_open, params->gap_extend, params->bandw,
					0, params->zdrop, s->len*params->match, &qle, &tle, &gtle, &gscore, &offset);

		if (gscore >= 0) { // to-end
			aln->read_start = 0;
			aln->ref_start = s->ref_pos - gtle;
			aln->truesc = gscore;
		} else {
			// TODO: did not reach the end of the query!
			return;
		}
	} else {
		aln->score = s->len * params->match;
		aln->read_start = 0;
		aln->ref_start = s->ref_pos;
	}

	if (s->read_pos + s->len != r->len) { // right extension
		int qle, tle, qe, re, gtle, gscore;
		int sc0 = aln->score;
		uint32 read_len = r->len - (s->read_pos + s->len);
		uint32 ref_len = hit_offset + hit_len - (s->ref_pos + s->len);
		int offset;
		aln->score = ksw_extend2(read_len, (const unsigned char*) r->seq.c_str(), ref_len, (const unsigned char*)ref.seq.c_str(), 5, params->score_matrix, params->gap_open, params->gap_extend, params->gap_open, params->gap_extend, params->bandw,
				0, params->zdrop, sc0, &qle, &tle, &gtle, &gscore, &offset);

		// similar to the above
		if(gscore >= 0) { // to-end extension
			aln->read_end = r->len;
			aln->ref_end = s->ref_pos + s->len + gtle;
			aln->truesc += gscore - sc0;
		} else {
			// TODO: did not reach the end of the query!
			return;
		}
	} else {
		aln->read_end = r->len;
		aln->ref_end = s->ref_pos + s->len;
	}
}

#define MAX_SEED_HITS 10

void process_read_hits_se(ref_t& ref, read_t* r, const index_params_t* params) {
	// index the read sequence: generate and store all kmers
	std::unordered_map<minhash_t, std::vector<uint32>> kmer2pos;
	for(uint32 i = 0; i <= (r->len - params->k); i++) {
		minhash_t kmer_hash = CityHash32(&r->seq[i], params->k);
		if(kmer2pos.find(kmer_hash) != kmer2pos.end()) {
			kmer2pos.emplace(std::make_pair(kmer_hash, std::vector<uint32>()));
		}
		kmer2pos[kmer_hash].push_back(i);
	}

	// 1. collect seeds shared between the reference and the read
	std::vector<seed_t> seeds;
	for(uint32 i = 0; i < r->ref_matches[r->best_n_hits].size(); i++) {
		ref_match_t ref_contig = r->ref_matches[r->best_n_hits][i];
		int hit_len = r->len + ref_contig.len;
		seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;
		uint32 step = 1;
		for(uint32 j = 0; j < hit_len; j += step) {
			minhash_t kmer_hash = CityHash32(&ref.seq[hit_offset + j], params->k);
			if(kmer2pos.find(kmer_hash) != kmer2pos.end()) { // found a shared kmer seed

				uint32 min_step = UINT_MAX;
				if(kmer2pos[kmer_hash].size() > MAX_SEED_HITS) continue; //skip if too frequent
				for(uint32 read_pos_idx = 0; read_pos_idx < kmer2pos[kmer_hash].size(); read_pos_idx++) {
					uint32 read_pos = kmer2pos[kmer_hash][read_pos_idx];
					bool true_hit = true;
					for(uint32 k = 0; k < params->k; k++) {
						if(ref.seq[hit_offset + j + k] != r->seq[read_pos + k]) {
							true_hit = false;
							break;
						}
					} if(!true_hit) continue;

					// extend as much as possible to the right
					uint32 e;
					for(e = params->k; e < hit_len; e++) {
						if(read_pos + e >= r->len) break;
						if(ref.seq[hit_offset + j + e] != r->seq[read_pos + e]) break;
					}
					seed_t s(hit_offset + j, read_pos, e, i);
					seeds.push_back(s);
					if(e < min_step) {
						min_step = e;
					}
				}
				// jump to next kmer position
				if(min_step != UINT_MAX) {
					step = min_step;
				} else {
					step = 1;
				}
			}
		}
	}
	if(seeds.size() == 0) return; // TODO: remove duplicates/collapse seeds

	// 2. chain seeds
	std::vector<chain_t> chains;
	kbtree_t(chn) *tree = kb_init(chn, KB_DEFAULT_SIZE);
	for (uint32 i = 0; i < seeds.size(); i++) {
		seed_t s = seeds[i];
		chain_t tmp;
		tmp.pos = s.ref_pos;
		bool new_chain = false;
		if (kb_size(tree)) {
			chain_t *lower, *upper;
			kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
			if (!lower || !test_and_merge(params, lower, &s)) new_chain = true;
		} else {
			new_chain = true;
		}
		if (new_chain) { // add the seed as a new chain
			tmp.seeds.push_back(s);
			kb_putp(chn, tree, &tmp);
		}
	}

	#define traverse_func(p_) (chains.push_back(*(p_)))
	__kb_traverse(chain_t, tree, traverse_func);
	#undef traverse_func

	// find the longest chain
	chain_t* max_chain = &chains[0];
	uint32 max_weight = get_chain_weight(max_chain);
	for (uint32 i = 1; i < chains.size(); i++) {
			chain_t* c = &chains[i];
			uint32 w = get_chain_weight(c);
			if(w > max_weight) {
				max_chain = c;
				max_weight = w;
			}
	}

	// 3. extend longest seeds with longest chains
	uint32 max_len = max_chain->seeds[0].len;
	seed_t* max_s = max_chain->seeds[0];
	for (uint32 i = 1; i < max_chain->seeds.size(); i++) {
		seed_t* s = &max_chain->seeds[i];
		if(s->len > max_len) {
			max_s = s;
			max_len = s->len;
		}
	}
	seed2alignment(max_s, ref, r, params);
	printf("Score: %u \n", r->aln.truesc);

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

			if(r->best_n_hits > 0) {
				process_read_hits_se(ref, r, params);
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
	//#pragma omp parallel for reduction(+:valid_hash, acc_hits, acc_top)
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		if(!reads.reads[i].valid_minhash) continue;
		//printf("Evaluating read %u \n", i);
		eval_read_hit(ref, &reads.reads[i], params);
		acc_hits += reads.reads[i].acc;
		acc_top += reads.reads[i].top_hit_acc;
		valid_hash += reads.reads[i].valid_minhash;
	}

	printf("Max number of windows matched by read %u \n", max_windows_matched);
	printf("Avg number of windows matched per read %.8f \n", (float) total_windows_matched/acc_hits);
	printf("Avg number of top contigs matched per read %.8f \n", (float) total_top_contigs/acc_hits);
	printf("Avg contig length per read %.8f \n", (float) total_contigs_length/total_windows_matched);
	printf("Avg diff of top 2 hits per read %.8f \n", (float) diff_num_top_hits/reads.reads.size());
	printf("Total number of accurate hits matching top = %d \n", acc_top);
	printf("Total number of accurate hits found = %d \n", acc_hits);
	printf("Total search time: %f sec\n", end_time - start_time);

}
