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
#include "sam.h"
#include "lsh.h"

static const int N_RAND_INLIERS = (getenv("N_RAND_INLIERS") ? atoi(getenv("N_RAND_INLIERS")) : 10);
static const int DELTA_POS = (getenv("DELTA_POS") ? atoi(getenv("DELTA_POS")) : 20);
static const bool MEDIAN = (getenv("MEDIAN") ? atoi(getenv("MEDIAN")) : false);
static const bool VERBOSE = (getenv("VERBOSE") ? atoi(getenv("VERBOSE")) : false);
static const int WEIGHT_INT = (getenv("WEIGHT_SCORE") ? atoi(getenv("WEIGHT_SCORE")) : 2);
static const bool WEIGHT_SCORE = (WEIGHT_INT == 1);
static const bool WEIGHT_SCORES_SEPARATELY = (WEIGHT_INT == 2);
static const bool WEIGHT_SCORES_MAX = (WEIGHT_INT == 3);
static const int CUTOFF = (getenv("CUTOFF") ? atoi(getenv("CUTOFF")) : 0);
static const float WEIGHT_SCALE = (getenv("WEIGHT_SCALE") ? atof(getenv("WEIGHT_SCALE")) : 1);


int compute_ref_contig_votes(ref_match_t ref_contig, ref_t& ref, read_t* r, const index_params_t* params);

struct heap_entry_t {
	seq_t pos;
	uint32_t len;
	uint16_t tid;
	uint16_t next_idx;
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

inline void heap_set_min(heap_entry_t* heap, int n) {
	int min_idx = 0;
	for(int i = 1; i < n; i++) {
		if(heap[i].pos < heap[min_idx].pos) {
			min_idx = i;
		}
	}
	heap_entry_t tmp = heap[0];
	heap[0] = heap[min_idx];
	heap[min_idx] = tmp;
}

inline void heap_update(heap_entry_t* heap, uint32 n) {
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

void sift_down(heap_entry_t* heap, int n, int i) {
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
	    sift_down(heap, n, min);
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
		heap[n-1] = tmp;
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

void process_merged_contig(seq_t contig_pos, uint32 contig_len, int n_diff_table_hits, ref_t& ref, read_t* r, const bool rc, const index_params_t* params) {
	// DEBUG
	if(r->ref_pos_l >= contig_pos - contig_len - params->ref_window_size && r->ref_pos_l <= contig_pos + params->ref_window_size) {
		r->collected_true_hit = true;
		r->processed_true_hit = true;
	}

	// filters
	if(contig_len > params->max_matched_contig_len) return;
	if(n_diff_table_hits < (int) params->min_n_hits) return;
	if(n_diff_table_hits < (int) (r->best_n_bucket_hits - params->dist_best_hit)) return;

	// passed filters
	if(n_diff_table_hits > r->best_n_bucket_hits) { // if more hits than best so far
		r->best_n_bucket_hits = n_diff_table_hits;
	}
	int contig_chunk_size = (r->len + 1000);
	int n_contigs = ceil(((double)contig_len/contig_chunk_size));
	for(int x = 0; x < n_contigs; x++) {
		seq_t p = contig_pos - (n_contigs-1-x)*contig_chunk_size;
		int l = contig_chunk_size;
		if(x == 0) {
			l = contig_len - (n_contigs-1)*contig_chunk_size;
		}
		ref_match_t rm(p, l, rc);
		compute_ref_contig_votes(rm, ref, r, params);
	}

	// DEBUG
	if(r->ref_pos_l >= contig_pos - contig_len - params->ref_window_size && r->ref_pos_l <= contig_pos + params->ref_window_size) {
		r->bucketed_true_hit = n_diff_table_hits;
	}

}

// output matches (ordered by the number of projections matched)
void collect_read_hits_contigs_inssort_pqueue(ref_t& ref, read_t* r, const bool rc, const index_params_t* params) {
	r->ref_matches.resize(params->n_tables);

	// priority heap of matched positions
	heap_entry_t heap[params->n_tables];
	int heap_size = 0;
	// push the first entries in each sorted bucket onto the heap
	for(uint32 t = 0; t < params->n_tables; t++) { // for each table
		if(r->ref_bucket_matches_by_table[t] == NULL) {
			continue;
		} else {
			heap[heap_size].pos = (*r->ref_bucket_matches_by_table[t])[0].pos;
			heap[heap_size].len = (*r->ref_bucket_matches_by_table[t])[0].len;
			heap[heap_size].tid = t;
			heap[heap_size].next_idx = 1;
		}
		heap_size++;
	}
	if(heap_size == 0) return; // all the matched buckets are empty
	heap_sort(heap, heap_size); // build heap

	int n_diff_table_hits = 0;
	int len = 0;
	seq_t last_pos = -1;
	std::bitset<256> occ;
	while(heap_size > 0) {
		heap_entry_t e = heap[0]; // get min
		seq_t e_last_pos = e.pos + e.len - 1;
		if(last_pos == (seq_t) -1 || (e.pos <= last_pos)) {
			if(last_pos == (seq_t) -1) { // first contig
				len = e.len - 1;
				last_pos = e_last_pos;
			} else if(last_pos < e_last_pos) { // extending contig
				len += e_last_pos - last_pos;
				last_pos = e_last_pos;
			}
			if(!occ.test(e.tid)) {
				n_diff_table_hits++;
			}
			occ.set(e.tid);
		} else {
			// found a boundary, store/handle last contig
			process_merged_contig(last_pos, len, n_diff_table_hits, ref, r, rc, params);

			// start a new contig
			n_diff_table_hits = 1;
			len = e.len - 1;
			last_pos = e_last_pos;
			occ.reset();
			occ.set(e.tid);
		}
		/*if(heap_size > 1) {
			seq_t next_min_pos_diff_bucket = heap[1].pos;
			if(e_last_pos < next_min_pos_diff_bucket) {
				while(e.next_idx < (*r->ref_bucket_matches_by_table[e.tid]).size()) {
					e_last_pos = (*r->ref_bucket_matches_by_table[e.tid])[e.next_idx].pos + (*r->ref_bucket_matches_by_table[e.tid])[e.next_idx].len - 1;
					if(e_last_pos < next_min_pos_diff_bucket) {
						e.next_idx++;
					} else {
						break;
					}
				}
			}
		} else {
			break; // only this bucket is left => remaining entries cannot have more than 1 hit
		}*/
		// push the next match from this bucket
		if(e.next_idx < (*r->ref_bucket_matches_by_table[e.tid]).size()) {
			heap[0].pos = (*r->ref_bucket_matches_by_table[e.tid])[e.next_idx].pos;
			heap[0].len = (*r->ref_bucket_matches_by_table[e.tid])[e.next_idx].len;
			heap[0].next_idx = e.next_idx+1;
			heap_update(heap, heap_size);
		} else { // no more entries in this bucket
			heap[0] = heap[heap_size - 1];
			heap_size--;
			heap_update(heap, heap_size);
		}
	}

	// add the last position
	if(last_pos != (seq_t) -1) {
		process_merged_contig(last_pos, len, n_diff_table_hits, ref, r, rc, params);
	}
}

// output matches (ordered by the number of projections matched)
void mark_contig_brackets(ref_t& ref, read_t* r, const bool rc, const index_params_t* params) {
	int bracket_size = params->ref_window_size/2;
	std::vector<int16_t>* ref_brackets;
	std::vector<int>* ref_brackets_dirty_mask;
	if(rc) {
		ref_brackets = r->ref_brackets_rc;
		ref_brackets_dirty_mask = r->ref_brackets_dirty_rc;
	} else {
		ref_brackets = r->ref_brackets_f;
		ref_brackets_dirty_mask = r->ref_brackets_dirty_f;
	}
	for(uint32 t = 0; t < params->n_tables; t++) { // for each table
		if(r->ref_bucket_matches_by_table[t] == NULL) {
			continue;
		} else {
			// for each window in the bucket
			for(uint32 i = 0; i < (*r->ref_bucket_matches_by_table[t]).size(); i++) {
				seq_t p = (*r->ref_bucket_matches_by_table[t])[i].pos;
				int len = (*r->ref_bucket_matches_by_table[t])[i].len;

				seq_t start_bracket = p/bracket_size;
				seq_t end_bracket = (p + params->ref_window_size + len - 1)/bracket_size;

				for(uint32 j = start_bracket; j < end_bracket; j++) {
					if((*ref_brackets_dirty_mask)[j] != (int) r->rid) {
						(*ref_brackets_dirty_mask)[j] = r->rid;
						(*ref_brackets)[j] = 0; // clear
					}
					if((*ref_brackets)[j] == -1) {
						continue; // this contig already was processed
					}
					(*ref_brackets)[j]++;

					int n_diff_table_hits = (*ref_brackets)[j];
					if(n_diff_table_hits >= (int) params->min_n_hits && n_diff_table_hits >= (int) (r->best_n_bucket_hits - params->dist_best_hit)) {
						if(n_diff_table_hits > r->best_n_bucket_hits) { // if more hits than best so far
							r->best_n_bucket_hits = n_diff_table_hits;
						}
						ref_match_t rm(j*bracket_size + bracket_size, bracket_size, rc);
						compute_ref_contig_votes(rm, ref, r, params);
						(*ref_brackets)[j] = -1; // mark as checked
					}
				}
			}
		}
	}
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
	uint32_t chain_end_pos = 0;
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
	ref_match_t ref_contig = r->ref_matches[r->best_n_bucket_hits][s.ref_contig_idx];
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

	for(uint32 i = 0; i < r->ref_matches[r->best_n_bucket_hits].size(); i++) {
		ref_match_t ref_contig = r->ref_matches[r->best_n_bucket_hits][i];
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
	chains.reserve(r->ref_matches[r->best_n_bucket_hits].size());
	for(uint32 i = 0; i < r->ref_matches[r->best_n_bucket_hits].size(); i++) {
		// REF CANDIDATE CONTIG
		ref_match_t ref_contig = r->ref_matches[r->best_n_bucket_hits][i];
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
	for(uint32 i = 0; i < r->ref_matches[r->best_n_bucket_hits].size(); i++) { // REF CANDIDATE CONTIG
		ref_match_t ref_contig = r->ref_matches[r->best_n_bucket_hits][i];
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


struct comp_shared_seeds
{
    bool operator()(const std::pair<minhash_t, uint32>& a, const std::pair<minhash_t, uint32>& b) const {
    	return a.second < b.second;
    }
};

int is_inform_kmer(const char* seq, const uint32_t len, const index_params_t* params) {
	uint32 base_counts[5] = { 0 };
	for(uint32 i = 0; i < len; i++) {
		base_counts[(int) seq[i]]++;
	}
	if(base_counts[4] > 0) { // N ambiguous bases
		return 0;
	}
	uint32 n_empty = 0;
	for(uint32 i = 0; i < 4; i++) {
		if(base_counts[i] == 0) {
			n_empty++;
		}
	}
	if(n_empty > 1) { // repetitions of 2 or 1 base
		return 0;
	}

	return 1;
}

// compute number of votes for this contig and its interior alignment
int compute_ref_contig_votes(ref_match_t ref_contig, ref_t& ref, read_t* r, const index_params_t* params) {
	const std::vector<std::pair<minhash_t, uint32>>& kmers = (ref_contig.rc) ? r->kmers_rc : r->kmers_f;

  //SDM pos == last position
	seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;
	seq_t padded_hit_offset = (hit_offset >= CONTIG_PADDING) ? hit_offset - CONTIG_PADDING : 0;
//SDM r == "read"
	uint32 search_len = ref_contig.len + 2*CONTIG_PADDING + r->len;
	std::vector<std::pair<minhash_t, uint32>> kmers_ref((search_len - params->k2 + 1));
	for(uint32 j = 0; j < search_len - params->k2 + 1; j++) {
			kmers_ref[j] = std::make_pair(ref.precomputed_kmer2_hashes[padded_hit_offset + j], padded_hit_offset+j);//std::make_pair(CityHash32(&ref.seq[padded_hit_offset + j], params->k2), padded_hit_offset+j);
	}
	std::sort(kmers_ref.begin(), kmers_ref.end());

	if(r->ref_pos_l >= ref_contig.pos - ref_contig.len - params->ref_window_size && r->ref_pos_l <= ref_contig.pos + params->ref_window_size) {
		//printf("RC %d contig pos %u len %u offset %u search_len %u \n", ref_contig.rc, ref_contig.pos, ref_contig.len, padded_hit_offset, search_len);
	}
	// find how many kmers are in common
	int kmer_votes = 0;
	uint64 aln_ref_pos = 0;
	const int n_rand_inliers = N_RAND_INLIERS;
	int rand_inliers_idx = 0;
	uint32 maybe_inliers[n_rand_inliers];
	uint64 avg_aln_pos = 0;
	const uint32 delta_pos = DELTA_POS;
	bool init_pass = true;
	bool median = MEDIAN;
	int votes_noransac = 0;
	int max_possible_votes_total = kmers.size();

	for(int p = 0; p < 2; p++) {
		uint32 idx_q = 0;
		uint32 idx_r = 0;
		while(idx_q < kmers.size() && idx_r < kmers_ref.size()) {
			uint32 kmer_hash_ref = kmers_ref[idx_r].first;
			uint32 kmer_hash_q = kmers[idx_q].first;
			if(kmer_hash_ref == kmer_hash_q) {
				// match
				uint32 match_aln_pos = kmers_ref[idx_r].second - kmers[idx_q].second;
				if(init_pass) {
					if(((idx_r < (kmers_ref.size()-1) && kmers_ref[idx_r + 1].first != kmer_hash_ref) || idx_r == kmers_ref.size()-1) &&
							((idx_r > 0 && kmers_ref[idx_r -1].first != kmer_hash_ref) || idx_r == 0) &&
							((idx_q < (kmers.size()-1) && kmers[idx_q + 1].first != kmer_hash_q) || idx_q == kmers.size()-1) &&
							((idx_q > 0 && kmers[idx_q -1].first != kmer_hash_q) || idx_q == 0)) { // unique kmer
						if(rand_inliers_idx < n_rand_inliers-1) {
							maybe_inliers[rand_inliers_idx] = match_aln_pos;
							rand_inliers_idx++;
						} else {
							maybe_inliers[rand_inliers_idx] = match_aln_pos;
							// find the average
							if (median)
							{
						        	std::sort(&maybe_inliers[0],maybe_inliers+n_rand_inliers);
								avg_aln_pos = maybe_inliers[(n_rand_inliers-1)/2];	

							}
							else
							{
								for(int z = 0; z < n_rand_inliers; z++) {
									avg_aln_pos += maybe_inliers[z];
								}
								avg_aln_pos = avg_aln_pos/n_rand_inliers;
							}

							init_pass = false;
							break;
						}
					}
				} else {
					// regardless, increment number of votes
                    votes_noransac++;
					if(match_aln_pos > (avg_aln_pos - delta_pos) && match_aln_pos < (avg_aln_pos + delta_pos)) {
						// within delta
						kmer_votes++;
						aln_ref_pos += match_aln_pos;
					}
				}
				idx_q++;
				idx_r++;
			} else if(kmers[idx_q].first < kmers_ref[idx_r].first) {
				idx_q++;
			} else {
				idx_r++;
			}
			int max_possible_votes = kmers.size() - idx_q + kmer_votes;
			if(max_possible_votes < r->max_votes_second_best) {// || max_possible_votes < 50) {
				break;  // if all the remaining votes cannot exceed max
			}
		}
	}

	// keep track of max and its alignment position
	if(kmer_votes > r->max_votes) {
		r->max_votes_second_best = r->max_votes;
		r->max_votes_noransac_second_best = r->max_votes_noransac;
		r->max_votes = kmer_votes;
		r->aln.ref_start = aln_ref_pos/r->max_votes;
		r->aln.rc = ref_contig.rc;
		r->max_votes_noransac = votes_noransac;
		r->max_possible_votes = max_possible_votes_total;

	} else if(kmer_votes > r->max_votes_second_best) {
		r->max_votes_second_best = kmer_votes;
		r->max_votes_noransac_second_best = votes_noransac;
	}

	if(r->ref_pos_l >= ref_contig.pos - ref_contig.len - params->ref_window_size && r->ref_pos_l <= ref_contig.pos + params->ref_window_size) {
		r->comp_votes_hit = kmer_votes;
	}
	return kmer_votes;
}

void process_read_hits_se_votes_opt(ref_t& ref, read_t* r, const index_params_t* params) {

	// index the read sequence: generate and store all kmers
	r->kmers_f.resize((r->len - params->k2 + 1));
	r->kmers_rc.resize((r->len - params->k2 + 1));
	for(uint32 i = 0; i < (r->len - params->k2 + 1); i++) {
		r->kmers_f[i] = std::make_pair(CityHash32(&r->seq[i], params->k2), i);
		r->kmers_rc[i] = std::make_pair(CityHash32(&r->rc[i], params->k2), i);
	}
	std::sort(r->kmers_f.begin(), r->kmers_f.end());
	std::sort(r->kmers_rc.begin(), r->kmers_rc.end());

	int n_top_buckets = params->n_top_buckets_search;
	int n_collected_hits = 0;
	int n_collected_buckets = 0;
	int idx = 0;
	while(n_collected_buckets < n_top_buckets) {
		if(r->best_n_bucket_hits - idx < 0) break;
		if(r->ref_matches[r->best_n_bucket_hits - idx].size() != 0) {
			n_collected_hits += r->ref_matches[r->best_n_bucket_hits - idx].size();
			n_collected_buckets++;
		}
		idx++;
	}
	std::vector<int> kmers_votes(n_collected_hits);
	std::vector<uint32> hit_bucket_index(n_collected_hits);
	std::vector<uint32> hit_bucket_pos(n_collected_hits);
	std::vector<uint64> aln_ref_pos(n_collected_hits);

	r->max_votes = 0;
	r->max_votes_second_best = 0;
	int top_contig_idx = 0;

	int n_proc_buckets = 0;
	idx = 0;
	int hit_idx = 0;
	while(n_proc_buckets < n_collected_buckets) {
		if(r->ref_matches[r->best_n_bucket_hits - idx].size() == 0) {
			idx++;
			continue;
		}
		for(uint32 i = 0; i < r->ref_matches[r->best_n_bucket_hits - idx].size(); i++) { // REF CANDIDATE CONTIG
			ref_match_t ref_contig = r->ref_matches[r->best_n_bucket_hits - idx][i];
			kmers_votes[hit_idx] = compute_ref_contig_votes(ref_contig, ref, r, params);

			// keep track of max
			if(kmers_votes[hit_idx] > r->max_votes) {
				r->max_votes_second_best = r->max_votes;
				r->max_votes = kmers_votes[hit_idx];
				top_contig_idx = hit_idx;
			} else if(kmers_votes[hit_idx] > r->max_votes_second_best) {
				r->max_votes_second_best = kmers_votes[hit_idx];
			}
			hit_bucket_index[hit_idx] = r->best_n_bucket_hits - idx;
			hit_bucket_pos[hit_idx] = i;
			hit_idx++;
		}
		idx++;
		n_proc_buckets++;
	}

	r->n_max_votes = 1 + (r->max_votes == r->max_votes_second_best ? 1 : 0);
	//if(r->max_votes != 0) {
	//	r->aln.ref_start = aln_ref_pos[top_contig_idx]/r->max_votes;
	//}


	if((r->max_votes > r->max_votes_second_best) && r->max_votes != 0) {// && max_votes > 50) {
		r->aln.score = 255*(r->max_votes - r->max_votes_second_best)/r->max_votes;
		printf("max_votes: %d, second_best: %d, score: %d\n", (int)r->max_votes, (int)r->max_votes_second_best, (int)r->aln.score);

		if(r->aln.score >= 30) {
			if(r->aln.score >= 30 && !(r->ref_pos_l >= r->aln.ref_start - 30 && r->ref_pos_l <= r->aln.ref_start + 30)) {
				if (VERBOSE)
				{
					printf("score %u max %u second %u true votes %u bucket %u max buckt %u true position  %u vs found %u \n", r->aln.score,
							r->max_votes, r->max_votes_second_best, r->comp_votes_hit, r->bucketed_true_hit, r->best_n_bucket_hits, r->ref_pos_l, r->aln.ref_start);
							}
			}
		}
	} else {
		r->aln.score = 0;
	}

	//seed_t s(first_kmer_match[top_contig_idx].first, first_kmer_match[top_contig_idx].second, params->k, top_contig_idx);
	//seed2alignment(s, ref, r, params);

	/*std::sort(kmers.begin(), kmers.end(), comp_shared_seeds()); // re-sort read kmers by position
	seq_t hit_offset = top_contig.pos - top_contig.len + 1;
	seq_t padded_hit_offset = (hit_offset >= CONTIG_PADDING) ? hit_offset - CONTIG_PADDING : 0;
	uint32 search_len = top_contig.len + 2*CONTIG_PADDING + r->len;

	int idx_q = 0;
	for(uint32 j = 0; j < search_len - params->k + 1; j++) {
		minhash_t h = CityHash32(&ref.seq[padded_hit_offset + j], params->k);
		while(idx_q < kmers.size() && (kmers[idx_q].first < h)) {
			idx_q++;
		}
		if(idx_q >= kmers.size()) break;
		if(kmers[idx_q].first == h) {
			// found match
			seed_t s(padded_hit_offset + j, kmers[idx_q].second, params->k, top_contig_idx);
			seed2alignment(s, ref, r, params);
			break;
		}
	}*/
}

struct comp_kmers
{
    bool operator()(const std::pair<minhash_t, uint32>& a, const std::pair<minhash_t, uint32>& b) const {
    	return a.first < b.first;
    }
};

void process_read_hits_se_votes_opt2(ref_t& ref, read_t* r, const index_params_t* params) {
	// index the read sequence: generate and store all kmers
	std::vector<std::pair<minhash_t, uint32>> kmers((r->len - params->k + 1));
	for(uint32 i = 0; i < (r->len - params->k + 1); i++) {
		kmers[i] = std::make_pair(CityHash32(&r->seq[i], params->k), i);
	}
	std::sort(kmers.begin(), kmers.end(), comp_kmers());

	std::vector<uint32> kmers_votes(r->ref_matches[r->best_n_bucket_hits].size());
	std::vector<std::pair<uint32, uint32>> first_kmer_match(r->ref_matches[r->best_n_bucket_hits].size());

	for(uint32 i = 0; i < r->ref_matches[r->best_n_bucket_hits].size(); i++) { // REF CANDIDATE CONTIG
		ref_match_t ref_contig = r->ref_matches[r->best_n_bucket_hits][i];
		seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;
		seq_t padded_hit_offset = (hit_offset >= CONTIG_PADDING) ? hit_offset - CONTIG_PADDING : 0;
		uint32 search_len = ref_contig.len + 2*CONTIG_PADDING + r->len;

		bool first_match = true;
		for(uint32 j = 0; j < search_len - params->k + 1; j++) {
			const minhash_t h = CityHash32(&ref.seq[padded_hit_offset + j], params->k);
			std::vector<std::pair<minhash_t, uint32>>::iterator it = std::lower_bound(kmers.begin(), kmers.end(), std::make_pair(h, (uint32) 0), comp_kmers());
			if(it != kmers.end() && !(h < (*it).first)) {
				kmers_votes[i]++; // found match
				if(first_match) {
					first_kmer_match[i] = std::make_pair(padded_hit_offset + j, (*it).second);
					first_match = false;
				}
			}
		}
	}
	uint32 top_contig_idx = std::distance(kmers_votes.begin(), std::max_element(kmers_votes.begin(), kmers_votes.end()));
	ref_match_t top_contig = r->ref_matches[r->best_n_bucket_hits][top_contig_idx];
	r->aln.ref_start = first_kmer_match[top_contig_idx].first - first_kmer_match[top_contig_idx].second;
}

#define MAX_TOP_HITS 100

void align_reads_minhash(ref_t& ref, reads_t& reads, const index_params_t* params) {
	printf("**** SRX Alignment: MinHash ****\n");

	for(uint32 x = 0; x < params->sketch_proj_indices.size(); x++) {
		printf("%u ", params->sketch_proj_indices[x]);
	}
	printf("\n");

	omp_set_num_threads(params->n_threads);
	//DEBUG--------
	#pragma omp parallel for
	for (uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		unsigned int pos_r;
		parse_read_mapping(r->name.c_str(), &r->seq_id, &r->ref_pos_l, &pos_r, &r->strand);
		r->seq_id = r->seq_id - 1;
		if(ref.subsequence_offsets.size() > 1) {
			r->ref_pos_l += ref.subsequence_offsets[r->seq_id]; // convert to global id
		}
	}//------------

	double start_time = omp_get_wtime();
	// split the reads across the threads
	#pragma omp parallel
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
			r->rid = i;

			// 1. index the read sequence: generate and store all k2 kmers
			r->kmers_f.resize((r->len - params->k2 + 1));
			r->kmers_rc.resize((r->len - params->k2 + 1));
			for(uint32 i = 0; i < (r->len - params->k2 + 1); i++) {
				r->kmers_f[i] = std::make_pair(CityHash32(&r->seq[i], params->k2), i);
				r->kmers_rc[i] = std::make_pair(CityHash32(&r->rc[i], params->k2), i);
			}
			std::sort(r->kmers_f.begin(), r->kmers_f.end());
			std::sort(r->kmers_rc.begin(), r->kmers_rc.end());

			// 2. index the read sequence using MinHash and collect the buckets
			if(!r->valid_minhash && !r->valid_minhash_rc) continue;
			r->ref_bucket_matches_by_table.resize(params->n_tables);
			if(r->valid_minhash) { // FORWARD
				for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
					minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes, params->sketch_proj_indices, t*params->sketch_proj_len);
					uint32 bucket_index = ref.hash_tables[t].bucket_indices[bucket_hash];
					if(bucket_index == ref.hash_tables[t].n_buckets) continue;
					if(ref.hash_tables[t].buckets_data_vectors[bucket_index].size() > 1000) {
						// DEBUG
						for(uint32 z = 0; z < ref.hash_tables[t].buckets_data_vectors[bucket_index].size(); z++) {
							seq_t p = ref.hash_tables[t].buckets_data_vectors[bucket_index][z].pos;
							uint32 len = ref.hash_tables[t].buckets_data_vectors[bucket_index][z].len;
							if(r->ref_pos_l >= p - len - params->ref_window_size && r->ref_pos_l <= p + params->ref_window_size) {
								r->collected_true_hit = true;
							}
						}
						continue;
					}
					r->any_bucket_hits = true;
					r->ref_bucket_matches_by_table[t] = &ref.hash_tables[t].buckets_data_vectors[bucket_index];
				}
				collect_read_hits_contigs_inssort_pqueue(ref, r, false, params);
			}
			if(r->valid_minhash_rc) { // RC
				for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
					minhash_t bucket_hash_rc = params->sketch_proj_hash_func.apply_vector(r->minhashes_rc, params->sketch_proj_indices, t*params->sketch_proj_len);
					uint32 bucket_index_rc = ref.hash_tables[t].bucket_indices[bucket_hash_rc];
					if(bucket_index_rc == ref.hash_tables[t].n_buckets) continue;
					if(ref.hash_tables[t].buckets_data_vectors[bucket_index_rc].size() > 1000) {
						// DEBUG
						for(uint32 z = 0; z < ref.hash_tables[t].buckets_data_vectors[bucket_index_rc].size(); z++) {
							seq_t p = ref.hash_tables[t].buckets_data_vectors[bucket_index_rc][z].pos;
							uint32 len = ref.hash_tables[t].buckets_data_vectors[bucket_index_rc][z].len;
							if(r->ref_pos_l >= p - len - params->ref_window_size && r->ref_pos_l <= p + params->ref_window_size) {
								r->collected_true_hit = true;
							}
						}
						continue; // no reference window fell into this bucket
					}
					r->any_bucket_hits = true;
					r->ref_bucket_matches_by_table[t] = &ref.hash_tables[t].buckets_data_vectors[bucket_index_rc];
				}
				collect_read_hits_contigs_inssort_pqueue(ref, r, true, params);
			}
			std::vector< VectorSeqPos* >().swap(r->ref_bucket_matches_by_table); //release memory

			if(!r->collected_true_hit) {
				VectorMinHash window_hashes(params->h);
				bool valid_hash = minhash(ref.seq.c_str(), r->ref_pos_l-1, params->ref_window_size,
							ref.high_freq_kmer_trie,
							ref.ignore_kmer_bitmask,
							marisa::Trie(), params,
							params->kmer_hasher, false,
							window_hashes);
				print_read(r);
				printf("READ HASHES F (valid = %d): \n", r->valid_minhash);
				for(uint32 x = 0; x < params->h; x++) {
					printf("%u ", r->minhashes[x]);
				}
				printf("\n");
				for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
					minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes, params->sketch_proj_indices, t*params->sketch_proj_len);
					printf("t = %u, hash = %u\n", t, bucket_hash);
					if(t == 0) {
						uint32 bucket_index = ref.hash_tables[t].bucket_indices[bucket_hash];
						for(uint32 z = 0; z < ref.hash_tables[t].buckets_data_vectors[bucket_index].size(); z++) {
							seq_t p = ref.hash_tables[t].buckets_data_vectors[bucket_index][z].pos;
							uint32 len = ref.hash_tables[t].buckets_data_vectors[bucket_index][z].len;
							printf("pos %u len %u \n", p, len);
						}	
					}
				}
				printf("\n");
				printf("READ HASHES RC (valid = %d): \n", r->valid_minhash_rc);
				for(uint32 x = 0; x < params->h; x++) {
					printf("%u ", r->minhashes_rc[x]);
				}
				printf("\n");
				for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
					minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes_rc, params->sketch_proj_indices, t*params->sketch_proj_len);
					printf("t = %u, hash = %u\n", t, bucket_hash);
				}
				printf("WINDOW: \n");
				for(uint32 x = r->ref_pos_l-1; x < r->ref_pos_l-1 + params->ref_window_size; x++) {
					printf("%c", iupacChar[(int)ref.seq[x]]);
				}
				printf("\n");
				printf("WINDOW HASHES (valid = %d): \n", valid_hash);
				for(uint32 x = 0; x < params->h; x++) {
					printf("%u ", window_hashes[x]);
				}
				printf("\n");
				for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
					minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(window_hashes, params->sketch_proj_indices, t*params->sketch_proj_len);
					printf("t = %u, hash = %u\n", t, bucket_hash);
				}

			}

			if(r->any_bucket_hits && (r->best_n_bucket_hits > 0)) {
				r->best_n_bucket_hits = r->best_n_bucket_hits - 1;
			}

			r->n_max_votes = 1 + (r->max_votes == r->max_votes_second_best ? 1 : 0);
			if((r->max_votes > r->max_votes_second_best) && r->max_votes > CUTOFF) {
				if (WEIGHT_SCORE) {
					r->aln.score = 250*(r->max_votes - r->max_votes_second_best)/(float)r->max_votes * r->max_votes_noransac/(float)r->max_possible_votes;
				}
				else if (WEIGHT_SCORES_SEPARATELY) {
					float second_best_score = 0;
					if(r->max_votes_noransac_second_best > 0) {	
						second_best_score =  r->max_votes_second_best*((float)r->max_votes_second_best/(float)r->max_votes_noransac_second_best);
					}
					r->aln.score = 250*(r->max_votes*((float)r->max_votes/(float)r->max_votes_noransac) - second_best_score)/(float)(r->max_votes*(float)r->max_votes/(float)r->max_votes_noransac);
				}
				else if (WEIGHT_SCORES_MAX) {
					r->aln.score = 30;
				}
				else {
					r->aln.score = 250*(r->max_votes - r->max_votes_second_best)/r->max_votes;
				}
				r->aln.score *= WEIGHT_SCALE;

				if (VERBOSE) {
					printf("max_votes: %d, second_best: %d, noransac: %d, possible: %d, score: %d\n", (int)r->max_votes, (int)r->max_votes_second_best, 
						(int)r->max_votes_noransac, (int)r->max_possible_votes, (int)r->aln.score);
					if(r->aln.score >= 30 && !(r->ref_pos_l >= r->aln.ref_start - 30 && r->ref_pos_l <= r->aln.ref_start + 30)) {
								printf("score %u max %u second %u true votes %u bucket %u max buckt %u true  %u found %u\n", r->aln.score,
										r->max_votes, r->max_votes_second_best, r->comp_votes_hit, r->bucketed_true_hit, r->best_n_bucket_hits, r->ref_pos_l, r->aln.ref_start);
										//,top_contig.pos, top_contig.len);
								//if(r->bucketed_true_hit) {
								//		print_read(r);
								//		printf("TRUE \n");
								//		for(uint32 x = r->ref_pos_l; x < r->ref_pos_l + 1000; x++) {
								//			printf("%c", iupacChar[(int)ref.seq[x]]);
								//		}
								//		printf("\n");
								//		printf("FALSE \n");
								//		for(uint32 x = r->aln.ref_start; x < r->aln.ref_start + 1000; x++) {
								//			printf("%c", iupacChar[(int)ref.seq[x]]);
								//		}
								//		printf("\n");
								//}


					}
				}
			} else {
				r->aln.score = 0;
			}
		}
	}

	// write the results to file
	store_alns_sam(reads, ref, params);

	double end_time = omp_get_wtime();
	printf("Collected read hits \n");

	printf("Evaluating read hits... \n");
	int valid_hash = 0;
	int mapped = 0;
	int n_collected = 0;
	int confident = 0;
	int acc_hits = 0;
	int acc_top = 0;
	int acc_dp = 0;
	int n_max_votes = 0;
	int best_hits = 0;
	int score = 0;
	int q10 = 0;
	int q30 = 0;
	int q30acc = 0;
	int q10acc = 0;
	int processed_true = 0;
	int bucketed_true = 0;
	int q30processed_true = 0;
	int q30bucketed_true = 0;
	//#pragma omp parallel for reduction(+:valid_hash, acc_hits, acc_top)
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->valid_minhash && !r->valid_minhash_rc) continue;
		valid_hash++;
		if(r->collected_true_hit) {
			n_collected++;
		}
		if(!r->any_bucket_hits) continue;
		if(r->processed_true_hit) {
			processed_true++;
		}
		if(r->bucketed_true_hit) {
			bucketed_true++;
		}
		if(r->aln.score <= 0){
			continue;
		}
		mapped++;
		//eval_read_hit(ref, r, params);
		acc_top += r->top_hit_acc;

		if(r->ref_pos_l >= r->aln.ref_start - 30 && r->ref_pos_l <= r->aln.ref_start + 30) {
			r->dp_hit_acc = 1;
		}
		acc_dp += r->dp_hit_acc;
		n_max_votes += r->n_max_votes;
		best_hits += r->best_n_bucket_hits + 1;
		score += r->aln.score;
		if(r->n_max_votes == 1) {
			confident++;
			if(r->dp_hit_acc) {
				acc_hits++;
			}
		}
		if(r->aln.score >= 30) {
			q30++;
			if(r->dp_hit_acc) {
				q30acc++;
			} else {
				//printf("score %u true  %u found %u \n", reads.reads[i].aln.score, reads.reads[i].ref_pos_l, reads.reads[i].aln.ref_start);
			}
			if(r->processed_true_hit) {
				q30processed_true++;
			}
			if(r->bucketed_true_hit) {
				q30bucketed_true++;
			}
		}
		if(r->aln.score >= 10) {
			q10++;
			if(r->dp_hit_acc) {
				q10acc++;
			}
		}
	}

	printf("Number of reads with valid F or RC hash %u \n", valid_hash);
	printf("Number of mapped reads score %u \n", mapped);
	printf("Number of mapped reads COLLECTED true hit %u \n", n_collected);
	printf("Number of mapped reads PROC true hit %u \n", processed_true);
	printf("Number of mapped reads BUCK true hit %u \n", bucketed_true);
	printf("Number of confidently mapped reads > 0 %u / accurate %u (%f pct)\n", confident, acc_hits, (float)acc_hits/(float)confident);
	printf("Number of confidently mapped reads Q10 %u / accurate %u (%f pct)\n", q10, q10acc, (float)q10acc/(float)q10);
	printf("Number of confidently mapped reads Q30 %u / accurate %u (%f pct)\n", q30, q30acc, (float)q30acc/(float)q30);
	printf("Number of confidently mapped reads Q30 PROCESSED true %u \n", q30processed_true);
	printf("Number of confidently mapped reads Q30 BUCKET true %u \n", q30bucketed_true);
	printf("Avg number of top votes matched per read %.8f \n", (float) n_max_votes/mapped);
	printf("Avg number of max bucket hits per read %.8f \n", (float) best_hits/mapped);
	printf("Avg score per read %.8f \n", (float) score/mapped);
	printf("Total search time: %f sec\n", end_time - start_time);

}
