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
#include <unordered_map>
#include <bitset>
#include <omp.h>
#include "index.h"
#include "io.h"
#include "hash.h"
#include "sam.h"
#include "lsh.h"

static const int N_INIT_ANCHORS_MAX = (getenv("N_INIT_ANCHORS_MAX") ? atoi(getenv("N_INIT_ANCHORS_MAX")) : 20);
//static const int N_TABLES_MAX = (getenv("N_TABLES_MAX") ? atoi(getenv("N_TABLES_MAX")) : 256);
static const int VERBOSE = (getenv("VERBOSE") ? atoi(getenv("VERBOSE")) : 0);

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

void generate_voting_kmers(const char* seq, const seq_t seq_offset, const seq_t seq_len,
		const uint8* kmer_mask, const uint16_t kmer_mask_len,
		const index_params_t* params,
		std::vector<std::pair<minhash_t, uint16_t>>& voting_kmers) {

	uint32 n_kmers = seq_len - kmer_mask_len + 1;
	voting_kmers.resize(n_kmers);
	for(uint16_t j = 0; j < n_kmers; j++) {
		/*std::string buffer;
		uint16_t n_set = 0;
		for(uint16_t k = 0; k < kmer_mask_len; k++) {
			if(kmer_mask[k]) {
				buffer.append(1, (int) seq[seq_offset + j + k]);
				n_set++;
			}
		}
		voting_kmers[j] = std::make_pair(CityHash32(buffer.c_str(), n_set), j);
		*/
		voting_kmers[j] = std::make_pair(CityHash32(&seq[seq_offset+j], kmer_mask_len), j);
		//std::make_pair(ref.precomputed_kmer2_hashes[padded_hit_offset + j], padded_hit_offset+j);
	}
}

inline bool pos_in_range_sig(int pos1, int pos2, uint32 delta) {
        return abs(pos1 - pos2) <= delta;
}


inline bool pos_in_range(uint32 pos1, uint32 pos2, uint32 delta) {
	return abs(pos1 - pos2) <= delta;
}

inline bool pos_in_range_asym(uint32 pos1, uint32 pos2, uint32 delta1, uint32 delta2) {
	uint32 delta11 = pos2 > delta1 ? delta1 : pos2;
	return pos1 > pos2 - delta11 && pos1 < pos2 + delta2;
}

inline bool is_unique_kmer(const std::vector<std::pair<minhash_t, uint16_t>>& sorted_kmers, const uint32 idx) {
	bool unique = true;
	minhash_t hash = sorted_kmers[idx].first;
	if(idx > 0 && sorted_kmers[idx-1].first == hash) {
		unique = false;
	} else if(idx < sorted_kmers.size()-1 && sorted_kmers[idx+1].first == hash) {
		unique = false;
	}
	return unique;
}


void votes_by_pos(std::vector<std::pair<minhash_t, uint16_t>> ref_kmers,
                std::vector<std::pair<minhash_t, uint16_t>> read_kmers, const index_params_t* params,
                uint32* n_inliers, int* pos, uint32* total_n_matches) {

	std::sort(ref_kmers.begin(), ref_kmers.end());
        std::sort(read_kmers.begin(), read_kmers.end());

	std::vector<int> votes(ref_kmers.size() + read_kmers.size());
	int pos0 = read_kmers.size();

	uint32 skip = 0;
	for(uint32 i = 0; i < read_kmers.size(); i++) {
		for(uint32 j = 0; j < ref_kmers.size(); j++) {
			// mark each repeat of the kmer
			if(read_kmers[i].first == ref_kmers[j].first) {
				int match_aln_pos = ref_kmers[j].second - read_kmers[i].second;
				//std::cout << pos0 + match_aln_pos << " " << pos0 << " " << match_aln_pos << " " << ref_kmers[i].second << " " << read_kmers[i].second << "\n"; 
				votes[pos0 + match_aln_pos]++;
			} //else if(ref_kmers[j].first > read_kmers[j].first) {
				//if(i < read_kmers.size()-1 && read_kmers[i+1].first > read_kmers[i].first) {
				//	skip = j;
				//}
			//	break;
			//}
		}

	}
	
	//std::cout << "Votes: \n";
	//for(int i = 0; i < votes.size(); i++) {
	//	std::cout << votes[i] << " ";
	//} std::cout << "\n";

	std::vector<int> votes_conv(ref_kmers.size() + read_kmers.size());
	int max = 0;
	uint32 max_idx = 0;
	for(uint32 i = 0; i < votes.size(); i++) {
		uint32 start = i > params->delta_inlier ? i - params->delta_inlier : 0;
		uint32 end = i + params->delta_inlier;
		if(end >= votes.size()) end = votes.size()-1;
		for(int j = start; j < end; j++) {
			votes_conv[i] += votes[j];
		}
		if(votes_conv[i] > max) {
			max = votes_conv[i];
			max_idx = i;
		}
	}

	//std::cout << "Votes Conv: \n";
        //for(int i = 0; i < votes.size(); i++) {
        //        std::cout << votes_conv[i] << " ";
        //} std::cout << "\n";
	//std::cout << "Max " << max << " max idx " << max_idx << "\n"; 

	if(max == 0) return;

	n_inliers[0] = max;
	pos[0] = max_idx - pos0;
	
	int max2 = 0;
	for(uint32 i = 0; i < votes_conv.size(); i++) {
		if(i == max_idx) continue;
		if(votes_conv[i] > max2 && !pos_in_range(max_idx, i, params->delta_x*params->delta_inlier)) {
			max2 = votes_conv[i];
			n_inliers[1] = max2;
			pos[1] = i - pos0;
			//break;
		}
	}
	//std::cout << "Max 2 " << n_inliers[1] << " max idx " << pos[1] << "\n";

}



#define UNIQUE1 0
#define UNIQUE2 1
void ransac(std::vector<std::pair<minhash_t, uint16_t>> ref_kmers,
		std::vector<std::pair<minhash_t, uint16_t>> read_kmers, const index_params_t* params,
		uint32* n_inliers, int* pos, uint32* total_n_matches) {
	std::sort(ref_kmers.begin(), ref_kmers.end());
	std::sort(read_kmers.begin(), read_kmers.end());

	int anchors[2][N_INIT_ANCHORS_MAX];
	uint32 n_anchors[2] = { 0 };
	int anchor_median_pos[2] = { 0 };

	for(int p = 0; p < 3; p++) {
		// start from the beginning in each pass
		uint32 idx_q = 0;
		uint32 idx_r = 0;
		while(idx_q < ref_kmers.size() && idx_r < read_kmers.size()) {
			if(ref_kmers[idx_r].first == read_kmers[idx_q].first) { // MATCH
				int match_aln_pos = ref_kmers[idx_r].second - read_kmers[idx_q].second;
				if(p == 0) { // -------- FIRST PASS: collect unique matches ---------
					if(is_unique_kmer(read_kmers, idx_r) && is_unique_kmer(ref_kmers, idx_q)) {
						if(n_anchors[UNIQUE1] < params->n_init_anchors) {
							anchors[UNIQUE1][n_anchors[UNIQUE1]] = match_aln_pos;
							n_anchors[UNIQUE1]++;
							//std::cout << " (1) " << ref_kmers[idx_r].second << " " << read_kmers[idx_q].second << "\n"; 
						} else { // found all the anchors
							break;
						}
					}
					idx_q++;
					idx_r++;
			} else if(p == 1) { // --------- SECOND PASS: collect DIFF unique matches -------
					if(n_anchors[UNIQUE1] == 0) break; // no unique matches exist
					// compute the median first pass position
					if(anchor_median_pos[UNIQUE1] == 0) {
						//std::sort(&anchors[UNIQUE1][0], &anchors[UNIQUE1][0] + n_anchors[UNIQUE1]);
						anchor_median_pos[UNIQUE1] = anchors[UNIQUE1][(n_anchors[UNIQUE1]-1)/2];
					}
					if(is_unique_kmer(read_kmers, idx_r) && is_unique_kmer(ref_kmers, idx_q)) {
						// if this position is not close to the first pick
						if(!pos_in_range_sig(match_aln_pos, anchor_median_pos[UNIQUE1], params->delta_x*params->delta_inlier)) {
							if(n_anchors[UNIQUE2] < params->n_init_anchors) {
								anchors[UNIQUE2][n_anchors[UNIQUE2]] = match_aln_pos;
								n_anchors[UNIQUE2]++;
								//std::cout << " (2) " << ref_kmers[idx_r].second << " " << read_kmers[idx_q].second << "\n";
							} else { // found all the anchors
								break;
							}
						}
					}
					idx_q++;
					idx_r++;
				} else { // --------- THIRD PASS: count inliers to each candidate position -------
					if(n_anchors[UNIQUE1] == 0) break;
					// find the median alignment positions
					if(n_anchors[UNIQUE2] > 0 && anchor_median_pos[UNIQUE2] == 0) {
						//std::sort(&anchors[UNIQUE2][0], &anchors[UNIQUE2][0] + n_anchors[UNIQUE2]);
						anchor_median_pos[UNIQUE2] = anchors[UNIQUE2][(n_anchors[UNIQUE2]-1)/2];
					}
					// count inliers
					for(int i = 0; i < 2; i++) {
						if(n_anchors[i] == 0) continue;
						if(pos_in_range_sig(anchor_median_pos[i], match_aln_pos, params->delta_inlier)) {
							n_inliers[i]++;
							pos[i] += match_aln_pos;
							//std::cout << " pos[" << i << "]=" << pos[i] << " match pos " << match_aln_pos << " " << ref_kmers[idx_r].second << " " <<  read_kmers[idx_q].second << "\n"; 
						}
					}
					if(idx_r < (ref_kmers.size()-1) && ref_kmers[idx_r + 1].first == ref_kmers[idx_r].first) {
						// if the next kmer is still a match in the reference, check it in next iteration as an inlier
						idx_r++;
					} else {
						(*total_n_matches)++;
						idx_q++;
						idx_r++;
					}
				}
			} else if(read_kmers[idx_q].first < ref_kmers[idx_r].first) {
				idx_q++;
			} else {
				idx_r++;
			}
		}
	}
	for(int i = 0; i < 2; i++) {
		if(n_inliers[i] > 0) {
			pos[i] = pos[i]/(int)n_inliers[i];
			//std::cout << " final pos " << pos[i] << " n_inliners " << n_inliers[i] << "\n";
		}
	}

	if(n_anchors[UNIQUE1] < 5) {
		n_inliers[0] = 0;
		n_inliers[1] = 0;
		votes_by_pos(ref_kmers, read_kmers, params, n_inliers, pos, total_n_matches);
	}
}


#define MATCH 1
#define ERROR 2
#define UNKNOWN 0

void propagate_matches(read_t* r, const int contig_pos,
		const std::vector<std::pair<minhash_t, uint16_t>>& contig_kmers,
		const std::vector<std::pair<minhash_t, uint16_t>>& read_kmers,
		const uint8* kmer_mask, const uint16_t kmer_mask_len,
		uint16_t* min_n_matches, uint16_t* min_n_errors,
		const index_params_t* params) {

	// find all the kmers in the read that have consistent matches
	std::vector<uint8> kmer_matches(read_kmers.size());
	uint16_t last_match_pos = contig_pos > (int) params->delta_inlier ? contig_pos - params->delta_inlier : 0;
	int found = 0;
	for(uint16_t i = 0; i < read_kmers.size(); i++) {
		minhash_t kmer = read_kmers[i].first;
		// find the match in the contig if exists
		int start_pos = contig_pos + i > (int) params->delta_inlier ? contig_pos + i - params->delta_inlier : 0;
		int end_pos = contig_pos + i + params->delta_inlier;
		if(end_pos < 0) end_pos = 0; 
		for(int j = start_pos; j < end_pos; j++) {
		//for(uint16_t j = 0; j < contig_kmers.size(); j ++) {
			if(j > (int) contig_kmers.size() - 1) break;
			if(contig_kmers[j].first == kmer) {
				kmer_matches[i] = 1;
				//last_match_pos = j; // the next match must occur at a later position
				found++;
				break;
			}
		}
	}

	//std::cout << found << "\n";

	// mark positions in the read that correspond to kmer matches
	std::vector<uint8> base_matches(r->len);
	for(uint16_t i = 0; i < kmer_matches.size(); i++) {
		if(kmer_matches[i] == 0) continue;
		for(uint16_t j = 0; j < kmer_mask_len; j++) {
			if(kmer_mask[j] == 1) {
				base_matches[i+j] = MATCH;
			}
		}
	}
	
	for(uint16_t i = 0; i < base_matches.size(); i++) {
		if(base_matches[i] == MATCH) {
			 (*min_n_matches) += 1;
		}
	}

	//std::cout << base_matches.size() << " matches " <<  *min_n_matches << "\n";
	// propagate known positions to determine errors and known positions
	for(uint16_t i = 0; i < kmer_matches.size(); i++) {
		if(kmer_matches[i]) continue;
		int n_kmer_bases = 0;
		int n_matched_kmer_bases = 0;
		uint16_t mismatch_pos = 0;
		for(uint16_t j = 0; j < kmer_mask_len; j++) {
			if(kmer_mask[j]) {
				n_kmer_bases++;
				if(base_matches[i+j] == MATCH) {
					n_matched_kmer_bases++;
				} else {
					mismatch_pos = i + j;
				}
			}
		}
		// must be an error if only one unknown
		if(n_matched_kmer_bases == n_kmer_bases - 1) {
			base_matches[mismatch_pos] = ERROR;
			(*min_n_errors)++;
		}
	}

	//std::cout << " inliers[0] " << kmer_inliers[0] << " inliers[1] " << kmer_inliers[1] << " min_match[0] " << min_match[0] << " min_match[1] " << min_match[1] << "\n";

}

#define KMER_MASK_LEN 20
//static const uint8 kmer_mask[KMER_MASK_LEN] = {1,0,1,1,0,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1};
//static const uint8 kmer_mask[KMER_MASK_LEN] =   {1,0,1,1,0,1,1,1,1,1,0,0,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,1};
static const uint8 kmer_mask[KMER_MASK_LEN] =     {1, 1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

// compute number of votes for this contig and its interior alignment
#define CONTIG_PADDING 100
int compute_ref_contig_votes(ref_match_t ref_contig, ref_t& ref, read_t* r, const index_params_t* params) {
	// get read voting kmers
	const std::vector<std::pair<minhash_t, uint16_t>>& read_kmers = (ref_contig.rc) ? r->kmers_rc : r->kmers_f;
	// generate ref voting kmers
	std::vector<std::pair<minhash_t, uint16_t>> ref_kmers;
	seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;
	seq_t padded_hit_offset = (hit_offset >= CONTIG_PADDING) ? hit_offset - CONTIG_PADDING : 0;
	uint32 search_len = ref_contig.len + 2*CONTIG_PADDING + r->len;
	generate_voting_kmers(ref.seq.c_str(), padded_hit_offset, search_len, kmer_mask, KMER_MASK_LEN, params, ref_kmers);	
	//std::sort(ref_kmers.begin(), ref_kmers.end());

	uint32 pos[2] = { 0 };
	int pos_sig[2] = { 0 };
	uint32 n_inliers[2] = { 0 };
	uint32 total_n_matches = 0;
	//ransac(ref_kmers, read_kmers, params, n_inliers, pos_sig, &total_n_matches);
	votes_by_pos(ref_kmers, read_kmers, params, n_inliers, pos_sig, &total_n_matches);

	//uint8 kmer_mask[20] = {1};
	for(int i = 0; i < 2; i++) {
		uint16_t min_n_matches = 0;
		uint16_t min_n_errors = 0;
		if(n_inliers[i] == 0) continue;
		propagate_matches(r, pos_sig[i], ref_kmers, read_kmers, kmer_mask, KMER_MASK_LEN, &min_n_matches, &min_n_errors, params);
		pos[i] = pos_sig[i] + padded_hit_offset;
		n_inliers[i] = 3*min_n_matches;
		//if(min_n_errors > 10) n_inliers[i] = 0;
	}

	// update votes and its alignment position
	for(int i = 0; i < 2; i++) {

		if(n_inliers[i] > r->top_aln.inlier_votes) {
			if(!pos_in_range(pos[i], r->top_aln.ref_start, 30)) {
				r->second_best_aln.inlier_votes = r->top_aln.inlier_votes;
				r->second_best_aln.total_votes = r->top_aln.total_votes;
				r->second_best_aln.ref_start = r->top_aln.ref_start;
			}
			// update best alignment
			r->top_aln.inlier_votes = n_inliers[i];
			r->top_aln.total_votes = total_n_matches;
			r->top_aln.ref_start = pos[i];
			r->top_aln.rc = ref_contig.rc;
		} else if(n_inliers[i] > r->second_best_aln.inlier_votes) {
			if(!pos_in_range(pos[i], r->top_aln.ref_start, 30)) {
				r->second_best_aln.inlier_votes = n_inliers[i];
				r->second_best_aln.total_votes = total_n_matches;
				r->second_best_aln.ref_start = pos[i];
			}
		}
	}
	//if(anchors_idx[UNIQUE1] < params->n_init_anchors) {
		// too few unique hits
	//	if(total_n_matches > r->max_total_votes_low_anchors) {
	//		r->max_total_votes_low_anchors = total_n_matches;
	//	}
	//}

	// keep track of max total votes
	if(total_n_matches > r->max_total_votes) {
		r->max_total_votes = total_n_matches;
	}

#if(SIM_EVAL)
	// DEBUG -----
	if(pos_in_range_asym(r->ref_pos_r, ref_contig.pos, ref_contig.len + params->ref_window_size, params->ref_window_size) ||
		pos_in_range_asym(r->ref_pos_l, ref_contig.pos, ref_contig.len + params->ref_window_size, params->ref_window_size)) {
		r->comp_votes_hit = n_inliers[0] > n_inliers[1] ? n_inliers[0] : n_inliers[1];
		//if(VERBOSE == 0) {
		//	printf("------True Hit------\n");
		//}
	}
	if(VERBOSE > 0) {
		printf("contig_rc = %d; contig_pos = %u; len = %u; aln_pos: %u %u; votes: %u %u; total matches: %u \n",
				ref_contig.rc, ref_contig.pos, ref_contig.len, pos[0], pos[1], n_inliers[0], n_inliers[1], total_n_matches);
	}
#endif
	return 0;
}

void process_merged_contig(seq_t contig_pos, uint32 contig_len, int n_diff_table_hits, ref_t& ref, read_t* r, const bool rc, const index_params_t* params) {
#if(SIM_EVAL)
	 if(pos_in_range_asym(r->ref_pos_r, contig_pos, contig_len + params->ref_window_size, params->ref_window_size) ||
		pos_in_range_asym(r->ref_pos_l, contig_pos, contig_len + params->ref_window_size, params->ref_window_size))  {
		r->collected_true_hit = true;
		r->processed_true_hit = true;
		r->true_n_bucket_hits = n_diff_table_hits;
	}
#endif

	// filters
	if(contig_len > params->max_matched_contig_len) return;
	if(n_diff_table_hits < (int) params->min_n_hits) return;
	//if(n_diff_table_hits < (int) (r->best_n_bucket_hits - params->dist_best_hit)) return;

	// passed filters
	if(n_diff_table_hits > r->best_n_bucket_hits) { // if more hits than best so far
		r->best_n_bucket_hits = n_diff_table_hits;
	}
	ref_match_t rm(contig_pos, contig_len, rc);
	compute_ref_contig_votes(rm, ref, r, params);

	/* ---- split the ref contig
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
	}*/

#if(SIM_EVAL)
	 if(pos_in_range_asym(r->ref_pos_r, contig_pos, contig_len + params->ref_window_size, params->ref_window_size) ||
		pos_in_range_asym(r->ref_pos_l, contig_pos, contig_len + params->ref_window_size, params->ref_window_size))  {
		r->bucketed_true_hit = n_diff_table_hits;
	}
#endif
}

int get_next_contig(const ref_t& ref, read_t* r, uint32 t, heap_entry_t* entry) {
	uint64 bid = r->ref_bucket_matches_by_table[t].first;
	if(bid == ref.index.bucket_offsets.size()) { // table ignored
		return 0;
	}

	minhash_t read_proj_hash = r->ref_bucket_matches_by_table[t].second;
	uint64 bucket_data_size = ref.index.bucket_offsets[bid+1] - ref.index.bucket_offsets[bid];
	uint64 bucket_data_offset = ref.index.bucket_offsets[bid];
	// get the next entry in the bucket that matches the read projection hash value
	while(entry->next_idx < bucket_data_size) {
		if(ref.index.buckets_data[bucket_data_offset + entry->next_idx].hash == read_proj_hash) {
			entry->pos = ref.index.buckets_data[bucket_data_offset + entry->next_idx].pos;
			entry->len = ref.index.buckets_data[bucket_data_offset + entry->next_idx].len;
			entry->tid = t;
			entry->next_idx++;
			return 1;
		}
		entry->next_idx++;
	}
	return 0;
}


#define N_TABLES_MAX 1024
// output matches (ordered by the number of projections matched)
void collect_read_hits(ref_t& ref, read_t* r, const bool rc, const index_params_t* params) {
	r->ref_matches.resize(params->n_tables);

	// priority heap of matched positions
	heap_entry_t heap[params->n_tables];
	int heap_size = 0;
	// push the first entries in each sorted bucket onto the heap
	for(uint32 t = 0; t < params->n_tables; t++) { // for each table
		heap[heap_size].next_idx = 0;
		if(get_next_contig(ref, r, t, &heap[heap_size]) > 0) {
			heap_size++;
		}
	}
	if(heap_size == 0) return; // all the matched buckets are empty
	heap_sort(heap, heap_size); // build heap

	int n_diff_table_hits = 0;
	int len = 0;
	seq_t last_pos = -1;
	std::bitset<N_TABLES_MAX> occ;
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
		// push the next match from this bucket
		if(get_next_contig(ref, r, e.tid, &heap[0]) > 0) {
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

#define MAX_BUCKET_SIZE 10000

void align_reads_minhash(ref_t& ref, reads_t& reads, const index_params_t* params) {
	printf("**** SRX Alignment: MinHash ****\n");
	omp_set_num_threads(params->n_threads);

#if(SIM_EVAL)
	#pragma omp parallel for
	for (uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		parse_read_mapping(r->name.c_str(), &r->seq_id, &r->ref_pos_l, &r->ref_pos_r, &r->strand);
		r->seq_id = r->seq_id - 1;
		if(ref.subsequence_offsets.size() > 1) {
			r->ref_pos_l += ref.subsequence_offsets[r->seq_id]; // convert to global id
			r->ref_pos_r += ref.subsequence_offsets[r->seq_id];
		}
	}
#endif

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

		int avg_score_per_thread = 0;
		int n_nonzero_scores = 0;
		for (uint32 i = chunk_start; i < chunk_end; i++) { // for each read of the thread's chunk
			if((i - chunk_start) % 1000 == 0 && (i - chunk_start) != 0) {
				printf("Thread %d processed %u reads \n", tid, i - chunk_start);
			}

			read_t* r = &reads.reads[i];
			r->rid = i;

			// 1. index the read sequence: generate and store all k2 kmers
			generate_voting_kmers(r->seq.c_str(), 0, r->len, kmer_mask, KMER_MASK_LEN, params, r->kmers_f);
			generate_voting_kmers(r->rc.c_str(), 0, r->len, kmer_mask, KMER_MASK_LEN, params, r->kmers_rc);
			//std::sort(r->kmers_f.begin(), r->kmers_f.end());
			//std::sort(r->kmers_rc.begin(), r->kmers_rc.end());

			// 2. index the read sequence using MinHash and collect the buckets
			if(!r->valid_minhash && !r->valid_minhash_rc) continue;
			r->ref_bucket_matches_by_table.resize(params->n_tables);
			if(r->valid_minhash) { // FORWARD
				for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
					minhash_t proj_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes, params->sketch_proj_indices, t*params->sketch_proj_len);
					minhash_t bucket_hash = params->sketch_proj_hash_func.bucket_hash(proj_hash);
					
					uint64_t bid = t*params->n_buckets + bucket_hash;
					uint32 bucket_size = ref.index.bucket_offsets[bid + 1] - ref.index.bucket_offsets[bid];
					uint64 bucket_data_offset = ref.index.bucket_offsets[bid];
					if(bucket_size >  MAX_BUCKET_SIZE) {
#if(SIM_EVAL)
						for(uint32 z = 0; z < bucket_size; z++) {
							seq_t p = ref.index.buckets_data[bucket_data_offset + z].pos;
							uint32 len = ref.index.buckets_data[bucket_data_offset + z].len;
							if(r->ref_pos_l >= p - len - params->ref_window_size && r->ref_pos_l <= p + params->ref_window_size) {
								r->collected_true_hit = true;
							}
						}
#endif
						r->ref_bucket_matches_by_table[t] = std::pair<uint64, minhash_t>(ref.index.bucket_offsets.size(), 0);
						continue;
					}
					r->any_bucket_hits = true;
					r->ref_bucket_matches_by_table[t] = std::pair<uint64, minhash_t>(bid, proj_hash);
				}
				collect_read_hits(ref, r, false, params);
			}
			if(r->valid_minhash_rc) { // RC
				for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
					minhash_t proj_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes_rc, params->sketch_proj_indices, t*params->sketch_proj_len);
					minhash_t bucket_hash = params->sketch_proj_hash_func.bucket_hash(proj_hash);

					uint64_t bid = t*params->n_buckets + bucket_hash;
					uint32 bucket_size = ref.index.bucket_offsets[bid + 1] - ref.index.bucket_offsets[bid];
					uint64 bucket_data_offset = ref.index.bucket_offsets[bid];
					if(bucket_size >  MAX_BUCKET_SIZE) {
#if(SIM_EVAL)
						for(uint32 z = 0; z < bucket_size; z++) {
							seq_t p = ref.index.buckets_data[bucket_data_offset + z].pos;
							uint32 len = ref.index.buckets_data[bucket_data_offset + z].len;
							if(r->ref_pos_l >= p - len - params->ref_window_size && r->ref_pos_l <= p + params->ref_window_size) {
								r->collected_true_hit = true;
							}
						}
#endif
						r->ref_bucket_matches_by_table[t] = std::pair<uint64, minhash_t>(ref.index.bucket_offsets.size(), 0);
						continue;
					}
					r->any_bucket_hits = true;
					r->ref_bucket_matches_by_table[t] = std::pair<uint64, minhash_t>(bid, proj_hash);
				}
				collect_read_hits(ref, r, true, params);
			}

			if(r->top_aln.inlier_votes > 0) {
				avg_score_per_thread += r->top_aln.inlier_votes;
				n_nonzero_scores++;
			}
		}
		if(n_nonzero_scores > 0) {	
			avg_score_per_thread = avg_score_per_thread/n_nonzero_scores;
		}
		// find the standard deviation
		/*uint32 std_dev = 0;
		for (uint32 i = chunk_start; i < chunk_end; i++) {
			read_t* r = &reads.reads[i];
			if(r->top_aln.inlier_votes > 0) {
				std_dev += pow((double)avg_score_per_thread - r->top_aln.inlier_votes, (double) 2);
			}
		}
		std_dev = sqrt((double)std_dev/n_nonzero_scores);
		int CUTOFF_3STDDEV = avg_score_per_thread - CUTOFF_STDDEV_WEIGHT*std_dev;
		if(CUTOFF_3STDDEV < 0) CUTOFF_3STDDEV = 0;*/

		// assign alignment quality scores
		for (uint32 i = chunk_start; i < chunk_end; i++) {
			read_t* r = &reads.reads[i];
			r->top_aln.score = 0;
			// top > 0 and top != second best
			if(r->top_aln.inlier_votes > r->second_best_aln.inlier_votes) {
				// if top contig has fewer non-inlier votes
				// only allow if there is a secondary hit, otherwise most likely an unanchored repeat
				//if(r->max_total_votes_low_anchors <= r->top_aln.total_votes) {

					// if sufficient votes were accumulated (lower thresholds for unique hit)
					if(r->top_aln.inlier_votes > params->votes_cutoff) {
						r->top_aln.score = 250*(r->top_aln.inlier_votes - r->second_best_aln.inlier_votes)/r->top_aln.inlier_votes;
						// scale by the distance from theoretical best possible votes
						if(avg_score_per_thread > 0 && params->enable_scale) {
							r->top_aln.score *= (float) r->top_aln.inlier_votes/(float) avg_score_per_thread;
						}

					}
				//}
				if(r->top_aln.rc) {
					r->top_aln.ref_start += r->len;
				}
			}

#if(SIM_EVAL)
			if(VERBOSE > 4 && !r->collected_true_hit) {
				print_read(r);
				printf("forward (valid = %d): \n", r->valid_minhash);
				for(uint32 x = 0; x < params->h; x++) {
					printf("%u ", r->minhashes[x]);
				}
				printf("\n");
				printf("rc (valid = %d): \n", r->valid_minhash_rc);
				for(uint32 x = 0; x < params->h; x++) {
					printf("%u ", r->minhashes_rc[x]);
				}
				printf("\n");
				printf("ref: \n");
				for(uint32 x = r->ref_pos_l-1; x < r->ref_pos_l-1 + params->ref_window_size; x++) {
					printf("%c", iupacChar[(int)ref.seq[x]]);
				}
				printf("\n");

			}
			if (VERBOSE ==  0 && r->top_aln.score >= 10 && !pos_in_range(r->ref_pos_r, r->top_aln.ref_start, 30) && !pos_in_range(r->ref_pos_l, r->top_aln.ref_start, 30)) {
				printf("WRONG: score %u max-votes: %u second-best-votes: %u true-contig-votes: %u true-bucket-hits: %u max-bucket-hits %u true-pos-l  %u true-pos-r: %u found-pos %u\n", 
						r->top_aln.score,
						r->top_aln.inlier_votes, r->second_best_aln.inlier_votes, r->comp_votes_hit, r->true_n_bucket_hits, r->best_n_bucket_hits,
						r->ref_pos_l, r->ref_pos_r, r->top_aln.ref_start);
				if(VERBOSE == 0 && r->true_n_bucket_hits == 0) { //r->bucketed_true_hit) {
					print_read(r);
					printf("true \n");
					for(uint32 x = r->ref_pos_l; x < r->ref_pos_l + r->len; x++) {
						printf("%c", iupacChar[(int)ref.seq[x]]);
					}
					printf("\n");
					seq_t p = r->top_aln.ref_start;
					if(r->top_aln.rc) {
						p -= r->len;
					}
					printf("false \n");
					for(uint32 x = p; x < p + r->len; x++) {
						printf("%c", iupacChar[(int)ref.seq[x]]);
					}
					printf("\n");
				}
			}
#endif
#if(DEBUG)
			if (VERBOSE > 0 && (pos_in_range(r->ref_pos_r, r->top_aln.ref_start, 30) || pos_in_range(r->ref_pos_l, r->top_aln.ref_start, 30))) {//ln.ref_start - 30 && r->ref_pos_l <= r->top_aln.ref_start + 30)) {
				printf("LOW: score %u max-votes: %u second-best-votes: %u true-contig-votes: %u true-bucket-hits: %u max-bucket-hits %u true-pos-l  %u true-pos-r: %u found-pos %u\n",
                                                r->top_aln.score,
                                                r->top_aln.inlier_votes, r->second_best_aln.inlier_votes, r->comp_votes_hit, r->true_n_bucket_hits, r->best_n_bucket_hits,
                                                r->ref_pos_l, r->ref_pos_r, r->top_aln.ref_start);
			}
#endif
		}
	}

	// write the results to file
	store_alns_sam(reads, ref, params);

	double end_time = omp_get_wtime();
	printf("Collected read hits \n");

#if(SIM_EVAL)
	printf("Evaluating read hits... \n");
	int valid_hash = 0;
	int n_collected = 0;
	int processed_true = 0;
	int bucketed_true = 0;
	int confident = 0;
	int acc = 0;
	int best_hits = 0;
	int true_pos_hits = 0;
	int score = 0;
	int max_votes_inl = 0;
	int max_votes_all = 0;
	int q10 = 0;
	int q30 = 0;
	int q30acc = 0;
	int q10acc = 0;
	int q30processed_true = 0;
	int q30bucketed_true = 0;
	int q10bucketed_true = 0;
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->valid_minhash && !r->valid_minhash_rc) continue;
		valid_hash++;
		
		true_pos_hits += r->true_n_bucket_hits;

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
		if(r->top_aln.score < 0){
			continue;
		}
		confident++;

		if(!r->top_aln.rc) {
			if(pos_in_range(r->ref_pos_l, r->top_aln.ref_start, 30)) {
				r->dp_hit_acc = 1;
			}
		} else {
			if(pos_in_range(r->ref_pos_r, r->top_aln.ref_start, 30)) {
				r->dp_hit_acc = 1;
			}
        }
		acc += r->dp_hit_acc;
		best_hits += r->best_n_bucket_hits;
		score += r->top_aln.score;
		max_votes_inl += r->top_aln.inlier_votes;
		max_votes_all += r->top_aln.total_votes;

		if(r->top_aln.score >= 30) {
			q30++;
			if(r->dp_hit_acc) {
				q30acc++;
			}
			if(r->processed_true_hit) {
				q30processed_true++;
			}
			if(r->bucketed_true_hit) {
				q30bucketed_true++;
			}
		}
		if(r->top_aln.score >= 10) {
			q10++;
			if(r->dp_hit_acc) {
				q10acc++;
			}
			if(r->bucketed_true_hit) {
				q10bucketed_true++;
			}
		}
	}
	float avg_true_pos_hits = (float) true_pos_hits/valid_hash;
	double stddev_true_pos_hits = 0;
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->valid_minhash && !r->valid_minhash_rc) continue;
		stddev_true_pos_hits += pow((double)avg_true_pos_hits - r->true_n_bucket_hits, (double) 2);
	}
	stddev_true_pos_hits = sqrt((double)stddev_true_pos_hits/valid_hash);

	printf("Number of reads with valid F or RC hash %u \n", valid_hash);
	printf("Number of mapped reads COLLECTED true hit %u \n", n_collected);
	printf("Number of mapped reads PROC true hit %u \n", processed_true);
	printf("Number of mapped reads BUCK true hit %u \n", bucketed_true);
	printf("Number of confidently mapped reads > 0 %u / accurate %u (%f pct)\n", confident, acc, (float)acc/(float)confident);
	printf("Number of confidently mapped reads Q10 %u / accurate %u (%f pct)\n", q10, q10acc, (float)q10acc/(float)q10);
	printf("Number of confidently mapped reads Q30 %u / accurate %u (%f pct)\n", q30, q30acc, (float)q30acc/(float)q30);
	printf("Number of confidently mapped reads Q30 PROCESSED true %u \n", q30processed_true);
	printf("Number of confidently mapped reads Q30 BUCKET true %u \n", q30bucketed_true);
	printf("Number of confidently mapped reads Q10 BUCKET true %u \n", q10bucketed_true);
	printf("Avg number of max bucket hits per read %.8f \n", (float) best_hits/confident);
	printf("MEAN number of bucket hits with TRUE pos per read %.8f \n", (float) avg_true_pos_hits);
	printf("STDDEV number of bucket hits with TRUE pos per read %.8f \n", (float) stddev_true_pos_hits);
	printf("Avg score per read %.8f \n", (float) score/confident);
	printf("Avg max inlier votes per read %.8f \n", (float) max_votes_inl/confident);
	printf("Avg max all votes per read %.8f \n", (float) max_votes_all/confident);
#endif
	printf("Total search time: %f sec\n", end_time - start_time);
}
