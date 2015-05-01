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
static const int N_TABLES_MAX = (getenv("N_TABLES_MAX") ? atoi(getenv("N_TABLES_MAX")) : 256);
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

// compute number of votes for this contig and its interior alignment
#define CONTIG_PADDING 100
int compute_ref_contig_votes(ref_match_t ref_contig, ref_t& ref, read_t* r, const index_params_t* params) {
	const std::vector<std::pair<minhash_t, uint32>>& kmers = (ref_contig.rc) ? r->kmers_rc : r->kmers_f;

	seq_t hit_offset = ref_contig.pos - ref_contig.len + 1;
	seq_t padded_hit_offset = (hit_offset >= CONTIG_PADDING) ? hit_offset - CONTIG_PADDING : 0;
	uint32 search_len = ref_contig.len + 2*CONTIG_PADDING + r->len;
	std::vector<std::pair<minhash_t, uint32>> kmers_ref((search_len - params->k2 + 1));
	for(uint32 j = 0; j < search_len - params->k2 + 1; j++) {
			kmers_ref[j] = std::make_pair(ref.precomputed_kmer2_hashes[padded_hit_offset + j], padded_hit_offset+j);//std::make_pair(CityHash32(&ref.seq[padded_hit_offset + j], params->k2), padded_hit_offset+j);
	}
	std::sort(kmers_ref.begin(), kmers_ref.end());

	// find how many kmers are in common
	int UNIQUE1 = 0;
	int UNIQUE2 = 1;
	int RAND = 2;
	uint32 sample_anchors[3][N_INIT_ANCHORS_MAX];
	int kmer_inliers[3] = { 0 };
	int total_kmer_matches = 0;
	uint32 anchors_idx[3] = { 0 };
	uint64 init_aln_pos[3] = { 0 };
	uint64 aln_ref_pos[3] = { 0 };

	for(int p = 0; p < 3; p++) {
		// start from the beginning in each pass
		uint32 idx_q = 0;
		uint32 idx_r = 0;
		while(idx_q < kmers.size() && idx_r < kmers_ref.size()) {
			uint32 kmer_hash_ref = kmers_ref[idx_r].first;
			uint32 kmer_hash_q = kmers[idx_q].first;
			if(kmer_hash_ref == kmer_hash_q) { // MATCH
				uint32 match_aln_pos = kmers_ref[idx_r].second - kmers[idx_q].second;

				// get the first (not necessarily unique) random matches
				if(anchors_idx[RAND] < params->n_init_anchors) {
					sample_anchors[RAND][anchors_idx[RAND]] = match_aln_pos;
					anchors_idx[RAND]++;
				}

				if(p == 0) { // -------- FIRST PASS: collect unique matches ---------
					if(((idx_r < (kmers_ref.size()-1) && kmers_ref[idx_r + 1].first != kmer_hash_ref) || idx_r == kmers_ref.size()-1) &&
							((idx_r > 0 && kmers_ref[idx_r -1].first != kmer_hash_ref) || idx_r == 0) &&
							((idx_q < (kmers.size()-1) && kmers[idx_q + 1].first != kmer_hash_q) || idx_q == kmers.size()-1) &&
							((idx_q > 0 && kmers[idx_q -1].first != kmer_hash_q) || idx_q == 0)) { // unique kmer
						if(anchors_idx[UNIQUE1] < params->n_init_anchors) {
							sample_anchors[UNIQUE1][anchors_idx[UNIQUE1]] = match_aln_pos;
							anchors_idx[UNIQUE1]++;
						} else { // found all the anchors
							break;
						}
					}
					idx_q++;
					idx_r++;
				} else if(p == 1) { // --------- SECOND PASS: collect DIFF unique matches -------
					if(anchors_idx[UNIQUE1] == 0) break; // no unique matches exist

					// compute the median first pass position
					if(init_aln_pos[UNIQUE1] == 0) {
						std::sort(&sample_anchors[UNIQUE1][0], &sample_anchors[UNIQUE1][0] + anchors_idx[UNIQUE1]);
						init_aln_pos[UNIQUE1] = sample_anchors[UNIQUE1][(anchors_idx[UNIQUE1]-1)/2];
					}
					if(((idx_r < (kmers_ref.size()-1) && kmers_ref[idx_r + 1].first != kmer_hash_ref) || idx_r == kmers_ref.size()-1) &&
							((idx_r > 0 && kmers_ref[idx_r -1].first != kmer_hash_ref) || idx_r == 0) &&
							((idx_q < (kmers.size()-1) && kmers[idx_q + 1].first != kmer_hash_q) || idx_q == kmers.size()-1) &&
							((idx_q > 0 && kmers[idx_q -1].first != kmer_hash_q) || idx_q == 0)) { // unique kmer
						// if this position is not close to the first pick
						if(match_aln_pos < init_aln_pos[UNIQUE1] - params->delta_x*params->delta_inlier || match_aln_pos > init_aln_pos[UNIQUE1] + params->delta_x*params->delta_inlier) {
							if(anchors_idx[UNIQUE2] < params->n_init_anchors) {
								sample_anchors[UNIQUE2][anchors_idx[UNIQUE2]] = match_aln_pos;
								anchors_idx[UNIQUE2]++;
							} else { // found all the anchors
								break;
							}
						}
					}
					idx_q++;
					idx_r++;
				} else { // --------- THIRD PASS: count inliers to each candidate position -------
					if(anchors_idx[UNIQUE1] == 0 && anchors_idx[RAND] == 0) break;

					// find the median alignment positions
					if(anchors_idx[UNIQUE2] > 0 && init_aln_pos[UNIQUE2] == 0) {
						std::sort(&sample_anchors[UNIQUE2][0], &sample_anchors[UNIQUE2][0] + anchors_idx[UNIQUE2]);
						init_aln_pos[UNIQUE2] = sample_anchors[UNIQUE2][(anchors_idx[UNIQUE2]-1)/2];
					}
					if(init_aln_pos[RAND] == 0) {
						std::sort(&sample_anchors[RAND][0], &sample_anchors[RAND][0] + anchors_idx[RAND]);
						init_aln_pos[RAND] = sample_anchors[RAND][(anchors_idx[RAND]-1)/2];
					}
					// count inliers
					for(int i = 0; i < 3; i++) {
						if(anchors_idx[i] == 0) continue;
						if(match_aln_pos > (init_aln_pos[i] - params->delta_inlier) && match_aln_pos < (init_aln_pos[i] + params->delta_inlier)) { // inliers (within delta)
							kmer_inliers[i]++;
							aln_ref_pos[i] += match_aln_pos;
						}
					}

					if(idx_r < (kmers_ref.size()-1) && kmers_ref[idx_r + 1].first == kmer_hash_ref) {
						// if the next kmer is still a match in the reference, check it in next iteration as an inlier
						idx_r++;
					} else {
						total_kmer_matches++;
						idx_q++;
						idx_r++;
					}
				}
			} else if(kmers[idx_q].first < kmers_ref[idx_r].first) {
				idx_q++;
			} else {
				idx_r++;
			}
			/*int max_possible_votes = kmers.size() - idx_q + kmer_votes;
			if(max_possible_votes < r->max_votes_second_best && max_possible_votes < r->max_votes_all) {
				break;  // if all the remaining votes cannot exceed max
			}*/
		}
	}

	for(int i = 0; i < 2; i++) {
		// keep track of max inlier votes and its alignment position
		if(kmer_inliers[i] > r->top_aln.inlier_votes) {
			if(aln_ref_pos[i]/kmer_inliers[i] < r->top_aln.ref_start - 30 || aln_ref_pos[i]/kmer_inliers[i] > r->top_aln.ref_start + 30) {
				r->second_best_aln.inlier_votes = r->top_aln.inlier_votes;
				r->second_best_aln.total_votes = r->top_aln.total_votes;
				r->second_best_aln.ref_start = r->top_aln.ref_start;
			}
			// update best alignment
			r->top_aln.inlier_votes = kmer_inliers[i];
			r->top_aln.total_votes = total_kmer_matches;
			r->top_aln.ref_start = aln_ref_pos[i]/kmer_inliers[i];
			r->top_aln.rc = ref_contig.rc;
		} else if(kmer_inliers[i] > r->second_best_aln.inlier_votes) {
			if(aln_ref_pos[i]/kmer_inliers[i] < r->top_aln.ref_start - 30 || aln_ref_pos[i]/kmer_inliers[i] > r->top_aln.ref_start + 30) {
				r->second_best_aln.inlier_votes = kmer_inliers[i];
				r->second_best_aln.total_votes = total_kmer_matches;
				r->second_best_aln.ref_start = aln_ref_pos[i]/kmer_inliers[i];
			}
		}
	}
	if(anchors_idx[UNIQUE1] < params->n_init_anchors) {
		// too few unique hits
		if(total_kmer_matches > r->max_total_votes_low_anchors) {
			r->max_total_votes_low_anchors = total_kmer_matches;
		}
	}

	// keep track of max total votes
	if(total_kmer_matches > r->max_total_votes) {
		r->max_total_votes = total_kmer_matches;
	}

#if(SIM_EVAL)
	// DEBUG -----
	if(r->ref_pos_l >= ref_contig.pos - ref_contig.len - params->ref_window_size && r->ref_pos_l <= ref_contig.pos + params->ref_window_size) {
		r->comp_votes_hit = kmer_inliers[0] > kmer_inliers[1] ? kmer_inliers[0] : kmer_inliers[1];
		if(VERBOSE > 1) {
			printf("RC %d contig pos %u len %u offset %u search_len %u \n", ref_contig.rc, ref_contig.pos, ref_contig.len, padded_hit_offset, search_len);
		}
	}
	if(VERBOSE > 3) {
		printf("aln pos %llu %llu %llu votes %u %u %u noransac votes %u avg pos %llu %llu %llu contig pos %u len %u \n",
						aln_ref_pos[0], aln_ref_pos[1], aln_ref_pos[2],
						kmer_inliers[0], kmer_inliers[1], kmer_inliers[2],
						total_kmer_matches, init_aln_pos[0], init_aln_pos[1], init_aln_pos[2], ref_contig.pos, ref_contig.len);
	}
#endif
	return 0;
}

void process_merged_contig(seq_t contig_pos, uint32 contig_len, int n_diff_table_hits, ref_t& ref, read_t* r, const bool rc, const index_params_t* params) {
#if(SIM_EVAL)
	if(r->ref_pos_l >= contig_pos - contig_len - params->ref_window_size && r->ref_pos_l <= contig_pos + params->ref_window_size) {
		r->collected_true_hit = true;
		r->processed_true_hit = true;
	}
#endif

	// filters
	if(contig_len > params->max_matched_contig_len) return;
	if(n_diff_table_hits < (int) params->min_n_hits) return;
	if(n_diff_table_hits < (int) (r->best_n_bucket_hits - params->dist_best_hit)) return;

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
	if(r->ref_pos_l >= contig_pos - contig_len - params->ref_window_size && r->ref_pos_l <= contig_pos + params->ref_window_size) {
		r->bucketed_true_hit = n_diff_table_hits;
	}
#endif
}

// output matches (ordered by the number of projections matched)
void collect_read_hits(ref_t& ref, read_t* r, const bool rc, const index_params_t* params) {
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

void align_reads_minhash(ref_t& ref, reads_t& reads, const index_params_t* params) {
	printf("**** SRX Alignment: MinHash ****\n");
	omp_set_num_threads(params->n_threads);

#if(SIM_EVAL)
	#pragma omp parallel for
	for (uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		unsigned int pos_r;
		parse_read_mapping(r->name.c_str(), &r->seq_id, &r->ref_pos_l, &pos_r, &r->strand);
		r->seq_id = r->seq_id - 1;
		if(ref.subsequence_offsets.size() > 1) {
			r->ref_pos_l += ref.subsequence_offsets[r->seq_id]; // convert to global id
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
#if(SIM_EVAL)
						for(uint32 z = 0; z < ref.hash_tables[t].buckets_data_vectors[bucket_index].size(); z++) {
							seq_t p = ref.hash_tables[t].buckets_data_vectors[bucket_index][z].pos;
							uint32 len = ref.hash_tables[t].buckets_data_vectors[bucket_index][z].len;
							if(r->ref_pos_l >= p - len - params->ref_window_size && r->ref_pos_l <= p + params->ref_window_size) {
								r->collected_true_hit = true;
							}
						}
#endif
						continue;
					}
					r->any_bucket_hits = true;
					r->ref_bucket_matches_by_table[t] = &ref.hash_tables[t].buckets_data_vectors[bucket_index];
				}
				collect_read_hits(ref, r, false, params);
			}
			if(r->valid_minhash_rc) { // RC
				for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
					minhash_t bucket_hash_rc = params->sketch_proj_hash_func.apply_vector(r->minhashes_rc, params->sketch_proj_indices, t*params->sketch_proj_len);
					uint32 bucket_index_rc = ref.hash_tables[t].bucket_indices[bucket_hash_rc];
					if(bucket_index_rc == ref.hash_tables[t].n_buckets) continue;
					if(ref.hash_tables[t].buckets_data_vectors[bucket_index_rc].size() > 1000) {
#if(SIM_EVAL)
						for(uint32 z = 0; z < ref.hash_tables[t].buckets_data_vectors[bucket_index_rc].size(); z++) {
							seq_t p = ref.hash_tables[t].buckets_data_vectors[bucket_index_rc][z].pos;
							uint32 len = ref.hash_tables[t].buckets_data_vectors[bucket_index_rc][z].len;
							if(r->ref_pos_l >= p - len - params->ref_window_size && r->ref_pos_l <= p + params->ref_window_size) {
								r->collected_true_hit = true;
							}
						}
#endif
						continue;
					}
					r->any_bucket_hits = true;
					r->ref_bucket_matches_by_table[t] = &ref.hash_tables[t].buckets_data_vectors[bucket_index_rc];
				}
				collect_read_hits(ref, r, true, params);
			}
			std::vector< VectorSeqPos* >().swap(r->ref_bucket_matches_by_table); //release memory

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
				if(r->max_total_votes_low_anchors <= r->top_aln.total_votes) {

					// if sufficient votes were accumulated (lower thresholds for unique hit)
					if(r->top_aln.inlier_votes > params->votes_cutoff) {

						r->top_aln.score = 250*(r->top_aln.inlier_votes - r->second_best_aln.inlier_votes)/r->top_aln.inlier_votes;

						// scale by the distance from theoretical best possible votes
						if(params->enable_scale) {
							r->top_aln.score *= (float) r->top_aln.inlier_votes/(float) avg_score_per_thread;
						}

					}
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
			if (VERBOSE > 0 && r->top_aln.score >= 30 && !(r->ref_pos_l >= r->top_aln.ref_start - 30 && r->ref_pos_l <= r->top_aln.ref_start + 30)) {
				printf("WRONG: score %u max %u second %u true votes %u bucket %u max buckt %u true  %u found %u\n", r->top_aln.score,
						r->top_aln.inlier_votes, r->second_best_aln.inlier_votes, r->comp_votes_hit, r->bucketed_true_hit, r->best_n_bucket_hits,
						r->ref_pos_l, r->top_aln.ref_start);
				if(VERBOSE > 3 && r->bucketed_true_hit) {
					print_read(r);
					printf("true \n");
					for(uint32 x = r->ref_pos_l; x < r->ref_pos_l + r->len; x++) {
						printf("%c", iupacChar[(int)ref.seq[x]]);
					}
					printf("\n");
					printf("false \n");
					for(uint32 x = r->top_aln.ref_start; x < r->top_aln.ref_start + r->len; x++) {
						printf("%c", iupacChar[(int)ref.seq[x]]);
					}
					printf("\n");
				}
			}
#endif
#if(DEBUG)
			if (VERBOSE > 0 && r->top_aln.score < 30 && r->top_aln.score >=10) { // && (r->ref_pos_l >= r->top_aln.ref_start - 30 && r->ref_pos_l <= r->top_aln.ref_start + 30)) {
					printf("LOW: score %u max %u second %u true votes %u bucket %u max buckt %u true  %u found %u second %u\n", r->top_aln.score,
									r->top_aln.inlier_votes, r->second_best_aln.inlier_votes, r->comp_votes_hit, r->bucketed_true_hit,
									r->best_n_bucket_hits, r->ref_pos_l, r->top_aln.ref_start, r->second_best_aln.ref_start);
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
	int score = 0;
	int max_votes_inl = 0;
	int max_votes_all = 0;
	int q10 = 0;
	int q30 = 0;
	int q30acc = 0;
	int q10acc = 0;
	int q30processed_true = 0;
	int q30bucketed_true = 0;
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
		if(r->top_aln.score <= 0){
			continue;
		}
		confident++;

		if(r->ref_pos_l >= r->top_aln.ref_start - 30 && r->ref_pos_l <= r->top_aln.ref_start + 30) {
			r->dp_hit_acc = 1;
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
		}
	}

	printf("Number of reads with valid F or RC hash %u \n", valid_hash);
	printf("Number of mapped reads COLLECTED true hit %u \n", n_collected);
	printf("Number of mapped reads PROC true hit %u \n", processed_true);
	printf("Number of mapped reads BUCK true hit %u \n", bucketed_true);
	printf("Number of confidently mapped reads > 0 %u / accurate %u (%f pct)\n", confident, acc, (float)acc/(float)confident);
	printf("Number of confidently mapped reads Q10 %u / accurate %u (%f pct)\n", q10, q10acc, (float)q10acc/(float)q10);
	printf("Number of confidently mapped reads Q30 %u / accurate %u (%f pct)\n", q30, q30acc, (float)q30acc/(float)q30);
	printf("Number of confidently mapped reads Q30 PROCESSED true %u \n", q30processed_true);
	printf("Number of confidently mapped reads Q30 BUCKET true %u \n", q30bucketed_true);
	printf("Avg number of max bucket hits per read %.8f \n", (float) best_hits/confident);
	printf("Avg score per read %.8f \n", (float) score/confident);
	printf("Avg max inlier votes per read %.8f \n", (float) max_votes_inl/confident);
	printf("Avg max all votes per read %.8f \n", (float) max_votes_all/confident);
#endif
	printf("Total search time: %f sec\n", end_time - start_time);
}
