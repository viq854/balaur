#include <stdlib.h>
#include <stdio.h>
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
#include "index.h"
#include "io.h"
#include "hash.h"
#include "cluster.h"

int find_window_match_diffk(ref_t& ref, cluster_t* cluster, const index_params_t* params);
int find_window_match(ref_t& ref, cluster_t* cluster, const index_params_t* params);
void find_windows_exact(ref_t& ref, read_t* r, const index_params_t* params);
void find_windows_exact_bucket(ref_t& ref, read_t* r, const index_params_t* params, int bucket);

int eval_cluster_hit(cluster_t* cluster);
int eval_read_hit(ref_t& ref, read_t* r);

void shuffle(int* perm);
void permute_ref(ref_t& ref, int perm[]);
void permute_reads(VectorClusters& reads, int perm[]);
void shift_bucket_ref(ref_t& ref, int bucket);
void shift_bucket_reads(reads_t& reads, int bucket);

void get_stats(const ref_t& ref, const VectorClusters& clusters) {
	uint32 printed = 0;
	for(uint32 i = 0; i < clusters.size(); i++) {
		cluster_t c = clusters[i];
		if(c.acc == 0) {
			if(printed > 50) break;
			printed++;
			printf("Cluster %u size = %zu hash %llx \n", i, c.reads.size(), c.simhash);
			printf("best diff = %d \n", c.best_hamd);
			//printf("best pos = %llu \n", c.ref_matches[c.best_pos]);
			for(uint32 j = 0; j < c.reads.size(); j++) {
				read_t* r = c.reads[j];
				print_read(r);
				unsigned int pos_l, pos_r;
				int strand;
				parse_read_mapping(r->name.c_str(), &pos_l, &pos_r, &strand);
//				unsigned int true_pos = pos_l - 1;
//				for(seq_t k = 0; k < ref.windows.size(); k++) {
//					if(true_pos == ref.windows[k].pos) {
//						printf("hamm dist true = %d simhash = %llx \n", hamming_dist(ref->windows[k].simhash, c.simhash), ref->windows[k].simhash);
//					}
//					if(c.ref_matches[c.best_pos] == ref->windows[k].pos) {
//						printf("hamm dist best = %d simhash = %llx\n", hamming_dist(ref->windows[k].simhash, c.simhash), ref->windows[k].simhash);
//					}
//				}
			}
		}
	}
}

// aligns the indexed reads to the indexed reference
// using multiple tables of read sampling
void align_reads_sampling(ref_t& ref, reads_t& reads, const index_params_t* params) {
	printf("**** SRX Alignment: Sampling ****\n");
	
//	// construct and search m hash tables using sampling
//	clock_t t = clock();
//	for(uint32 i = 0; i < params->m; i++) {
//		printf("Table %d \n", i);
//		index_ref_table_i(ref, params, i);
//		index_reads_table_i(reads, params, i);
//
//		for(uint32 j = 0; j < reads.reads.size(); j++) {
//			if(reads.reads[j].acc == 1) continue;
//			// binary search to find the matching ref window(s)
//			find_windows_exact(ref, &reads.reads[j], params);
//			eval_read_hit(&reads.reads[j]);
//		}
//	}
//
//	int acc_hits = 0;
//	for(uint32 i = 0; i < reads.reads.size(); i++) {
//		acc_hits += reads.reads[i].acc;
//	}
//	printf("Total number of accurate hits found = %d \n", acc_hits);
//	printf("Total search time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
}

int eval_minhash_hits_old(read_t* r, const index_params_t* params) {
//   unsigned int pos_l, pos_r;
//   int strand;
//   parse_read_mapping(r->name.c_str(), &pos_l, &pos_r, &strand);
//
//   for(seq_t j = pos_l - 10; j <= pos_r + 10; j++) {
//	   MapPos2MinCount::const_iterator v;
//	   if((v = r->matched_window_counts.find(j)) != r->matched_window_counts.end()) {
//		  if(v->second > params->n_min_matched) {
//			  r->acc = 1;
//			  break;
//		  }
//	   }
//
//    	if(r->acc == 1) {
//    		break;
//    	}
//    }
    return (r->acc == 1);
}

//int eval_read_hit(ref_t& ref, read_t* r) {
//   unsigned int pos_l, pos_r;
//   int strand;
//   parse_read_mapping(r->name.c_str(), &pos_l, &pos_r, &strand);
//
//   for(seq_t j = pos_l - 20; j <= pos_r + 20; j++) {
//	   for(uint32 t = 0; t < r->ref_bucket_id_matches.size(); t++) {
//		   buckets_t* buckets = &ref.hash_tables[t];
//		   for(seq_t idx = 0; idx < r->ref_bucket_id_matches[t].size(); idx++) {
//			   uint32 bucket_index = r->ref_bucket_id_matches[t][idx];
//			   VectorSeqPos& bucket = buckets->buckets_data_vectors[bucket_index];
//			   for(uint32 match = 0; match < bucket.size(); match++) {
//				   seq_t hit_pos = bucket[match];
//				   if(hit_pos == j) {
//					   r->acc = 1;
//					   break;
//				   }
//			   }
//			   if(r->acc == 1) {
//				   break;
//			   }
//		   }
//		   if(r->acc == 1) {
//		       	break;
//		   }
//        }
//    	if(r->acc == 1) {
//    		break;
//    	}
//    }
//    return (r->acc == 1);
//}

int eval_read_hit(ref_t& ref, read_t* r, const index_params_t* params) {
   unsigned int pos_l, pos_r;
   int strand;
   parse_read_mapping(r->name.c_str(), &pos_l, &pos_r, &strand);

   r->top_hit_acc = 0;
   
   //printf("best_n %u \n", r->best_n_hits);
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


void process_read_hits(ref_t& ref, read_t* r, const index_params_t* params) {

	// find the contig with most hits

	// find the contig with second most hits
}

void collect_read_hits_all(ref_t& ref, read_t* r, const index_params_t* params) {

	// construct a priority heap of matches ordered by the number of projections matched
	r->ref_matches.resize(params->n_tables);

	uint32 n_best_hits = 0; // best number of table hits found so far
	VectorU32 bucket_data_consumed_indices(params->n_tables);

	// process all the hits in intervals
	std::vector<seq_t> matches;
	matches.reserve(1000);
	for(uint32 i = 1; i <= ceil((float) ref.len/params->hit_collection_interval); i++) {
		matches.clear();
		for(uint32 t = 0; t < params->n_tables; t++) { // for each table
			buckets_t* buckets = &ref.hash_tables[t];
			uint32 bucket_index = r->ref_bucket_id_matches_by_table[t];
			if(bucket_index == buckets->n_buckets) {
				continue; // no reference window fell into this bucket
			}
			VectorSeqPos& bucket = buckets->buckets_data_vectors[bucket_index];
			uint32 data_pointer = bucket_data_consumed_indices[t];
			for(uint32 match = data_pointer; match < bucket.size(); match++) {
				if(bucket[match] < i*params->hit_collection_interval) {
					matches.push_back(bucket[match]);
					bucket_data_consumed_indices[t]++;
				} else {
					break; // requires that each bucket is sorted
				}
			}
		}
		
		if(matches.size() == 0) continue;

		// count how many times a position occurs
		std::sort(matches.begin(), matches.end());

		seq_t last_pos = matches[0];
		uint32 n_diff_table_hits = 1;
		for(uint32 i = 1; i < matches.size(); i++) {
			seq_t pos = matches[i];
			if(pos == last_pos) {
				n_diff_table_hits++;
			} else {
				// found a boundary, store
				if(n_diff_table_hits >= params->min_n_hits) {
					if(n_diff_table_hits > n_best_hits) { // if more hits than best so far
						n_best_hits = n_diff_table_hits;
						if(r->ref_matches[n_diff_table_hits-1].size() < params->max_best_hits) {
							ref_match_t rm(last_pos, 0);
							r->ref_matches[n_diff_table_hits-1].push_back(rm);
						}
					} else {
						// sub-optimal
						// only store if the number of hits is not too much lower than the best so far
						if(n_best_hits < params->dist_best_hit || n_diff_table_hits > (n_best_hits - params->dist_best_hit)) {
							if(r->ref_matches[n_diff_table_hits-1].size() < params->max_suboptimal_hits) {
								ref_match_t rm(last_pos, 0);
								r->ref_matches[n_diff_table_hits-1].push_back(rm);
							}
						}
					}
				}
				n_diff_table_hits = 1;
			}
			last_pos = pos;
		}

		// add the last position
		if(n_diff_table_hits >= params->min_n_hits) {
			if(n_diff_table_hits > n_best_hits) { // if more hits than best so far
				n_best_hits = n_diff_table_hits;
				if(r->ref_matches[n_diff_table_hits-1].size() < params->max_best_hits) {
					ref_match_t rm(last_pos, 0);
					r->ref_matches[n_diff_table_hits-1].push_back(rm);
				}
			} else {
				// sub-optimal
				// only store if the number of hits is not too much lower than the best so far
				if(n_best_hits < params->dist_best_hit || n_diff_table_hits > (n_best_hits - params->dist_best_hit)) {
					if(r->ref_matches[n_diff_table_hits-1].size() < params->max_suboptimal_hits) {
						ref_match_t rm(last_pos, 0);
						r->ref_matches[n_diff_table_hits-1].push_back(rm);
					}
				}
			}
		}

	}
	r->best_n_hits = (n_best_hits > 0) ? n_best_hits - 1 : 0;
	r->ref_bucket_id_matches_by_table = VectorU32(); //release memory
}

void collect_read_hits_contigs(ref_t& ref, read_t* r, const index_params_t* params) {

	// construct a priority heap of matches ordered by the number of projections matched
	r->ref_matches.resize(params->n_tables);

	uint32 n_best_hits = 0; // best number of table hits found so far
	VectorU32 bucket_data_consumed_indices(params->n_tables);

	// process all the hits in intervals
	std::vector<std::pair<seq_t, uint32> > matches;
	matches.reserve(1000);
	for(uint32 i = 1; i <= ceil((float) ref.len/params->hit_collection_interval); i++) {
		matches.clear();
		for(uint32 t = 0; t < params->n_tables; t++) { // for each table
			buckets_t* buckets = &ref.hash_tables[t];
			uint32 bucket_index = r->ref_bucket_id_matches_by_table[t];
			if(bucket_index == buckets->n_buckets) {
				continue; // no reference window fell into this bucket
			}
			VectorSeqPos& bucket = buckets->buckets_data_vectors[bucket_index];
			uint32 data_pointer = bucket_data_consumed_indices[t];
			for(uint32 match = data_pointer; match < bucket.size(); match++) {
				if(bucket[match] < i*params->hit_collection_interval) {
					matches.push_back(std::make_pair(bucket[match], t));
					bucket_data_consumed_indices[t]++;
				} else {
					break; // requires that each bucket is sorted
				}
			}
		}
		if(matches.size() == 0) continue;

		// count how many times a position occurs
		std::sort(matches.begin(), matches.end());

		seq_t last_pos = matches[0].first;
		VectorBool occ(params->n_tables, false);
		occ[matches[0].second] = true;
		uint32 n_diff_table_hits = 1;
		uint32 len = 0;
		for(uint32 i = 1; i < matches.size(); i++) {
			seq_t pos = matches[i].first;
			if(pos <= last_pos + params->contig_gap) {
				if(!occ[matches[i].second]) {
					n_diff_table_hits++;
				}
				occ[matches[i].second] = true;
				len += pos - last_pos;
			} else {
				// found a boundary, store
				assert(n_diff_table_hits > 0);
				assert(n_diff_table_hits-1 < params->n_tables);
				if(n_diff_table_hits >= params->min_n_hits) {
					if(n_diff_table_hits > n_best_hits) { // if more hits than best so far
						n_best_hits = n_diff_table_hits;
						if(r->ref_matches[n_diff_table_hits-1].size() < params->max_best_hits) {
							ref_match_t rm(last_pos, len);
							r->ref_matches[n_diff_table_hits-1].push_back(rm);
						}
					} else {
						// sub-optimal
						// only store if the number of hits is not too much lower than the best so far
						if(n_best_hits < params->dist_best_hit || n_diff_table_hits > (n_best_hits - params->dist_best_hit)) {
							if(r->ref_matches[n_diff_table_hits-1].size() < params->max_suboptimal_hits) {
								ref_match_t rm(last_pos, len);
								r->ref_matches[n_diff_table_hits-1].push_back(rm);
							}
						}
					}
				}
				std::fill(occ.begin(), occ.end(), false);
				occ[matches[i].second] = true;
				n_diff_table_hits = 1;
				len = 0;
			}
			last_pos = pos;
		}

		// add the last position
		assert(n_diff_table_hits > 0);
		assert(n_diff_table_hits-1 < params->n_tables);
		if(n_diff_table_hits >= params->min_n_hits) {
			if(n_diff_table_hits > n_best_hits) { // if more hits than best so far
				n_best_hits = n_diff_table_hits;
				if(r->ref_matches[n_diff_table_hits-1].size() < params->max_best_hits) {
					ref_match_t rm(last_pos, len);
					r->ref_matches[n_diff_table_hits-1].push_back(rm);
				}
			} else {
				// sub-optimal
				// only store if the number of hits is not too much lower than the best so far
				if(n_best_hits < params->dist_best_hit || n_diff_table_hits > (n_best_hits - params->dist_best_hit)) {
					if(r->ref_matches[n_diff_table_hits-1].size() < params->max_suboptimal_hits) {
						ref_match_t rm(last_pos, len);
						r->ref_matches[n_diff_table_hits-1].push_back(rm);
					}
				}
			}
		}

	}
	r->best_n_hits = (n_best_hits > 0) ? n_best_hits - 1 : 0;
	assert(r->best_n_hits < params->n_tables);
	r->ref_bucket_id_matches_by_table = VectorU32(); //release memory
}

//int comp_clusters(const cluster_t r1, const cluster_t r2) {
//	hash_t h1 = r1.simhash;
//	hash_t h2 = r2.simhash;
//	return (h1 > h2) - (h1 < h2);
//}

void collect_read_hits_contigs_inplace_merge(ref_t& ref, read_t* r, const index_params_t* params) {

	// output matches (ordered by the number of projections matched)
	r->ref_matches.resize(params->n_tables);
	uint32 n_best_hits = 0; // best number of table hits found so far

	// construct a priority heap of matched positions
	// comparing by position first
	typedef std::pair<seq_t, uint32> PairPosTid;
	std::priority_queue<PairPosTid, std::vector<PairPosTid>, std::greater<PairPosTid> > q;
	VectorU32 bucket_indices(params->n_tables);

	// push the first entries in each sorted bucket onto the heap
	for(uint32 t = 0; t < params->n_tables; t++) { // for each table
		buckets_t* buckets = &ref.hash_tables[t];
		uint32 bucket_index = r->ref_bucket_id_matches_by_table[t];
		printf("Table %u - index %u \n", t, bucket_index);
		if(bucket_index == buckets->n_buckets) {
			continue; // no reference window fell into this bucket
		}
		VectorSeqPos& bucket = buckets->buckets_data_vectors[bucket_index];
		q.push(std::make_pair(bucket[0], t));
		printf("Pushed match (%u, %u) \n", bucket[0], t);
		bucket_indices[t]++;
	}
	if(q.empty()) return; // all the matched buckets are empty

	VectorBool occ(params->n_tables, false);
	PairPosTid top = q.top();
	q.pop();
	uint32 tid = top.second;
	// push the next match from this bucket
	VectorSeqPos& bucket = ref.hash_tables[tid].buckets_data_vectors[r->ref_bucket_id_matches_by_table[tid]];
	if(bucket_indices[tid] < bucket.size()) {
		q.push(std::make_pair(bucket[bucket_indices[tid]], tid));
		printf("Pushed match (%u, %u) \n", bucket[bucket_indices[tid]], tid);
		bucket_indices[tid]++;
	}

	seq_t last_pos = top.first;
	occ[tid] = true;
	uint32 n_diff_table_hits = 1;
	uint32 len = 0;
	while(!q.empty()) {
		top = q.top();
		q.pop();
		seq_t pos = top.first;
		uint32 tid = top.second;
		if(pos <= last_pos + params->contig_gap) {
			if(!occ[tid]) {
				n_diff_table_hits++;
			}
			occ[tid] = true;
			len += pos - last_pos;
			printf("Pop within contig cur_pos %u last_pos %u tid %u len %u\n", pos, last_pos, tid, len);
		} else { // found a boundary, store
			printf("Pop found boundary cur_pos %u last_pos %u tid %u len %u\n", pos, last_pos, tid, len);
			assert(n_diff_table_hits > 0);
			assert(n_diff_table_hits-1 < params->n_tables);
			if(n_diff_table_hits >= params->min_n_hits) {
				if(n_diff_table_hits > n_best_hits) { // if more hits than best so far
					n_best_hits = n_diff_table_hits;
					if(r->ref_matches[n_diff_table_hits-1].size() < params->max_best_hits) {
						ref_match_t rm(last_pos, len);
						r->ref_matches[n_diff_table_hits-1].push_back(rm);
					}
				} else {
					// sub-optimal
					// only store if the number of hits is not too much lower than the best so far
					if(n_best_hits < params->dist_best_hit || n_diff_table_hits > (n_best_hits - params->dist_best_hit)) {
						if(r->ref_matches[n_diff_table_hits-1].size() < params->max_suboptimal_hits) {
							ref_match_t rm(last_pos, len);
							r->ref_matches[n_diff_table_hits-1].push_back(rm);
						}
					}
				}
				printf("MATCH pos %u len %u n_hits %u \n", pos, len, n_diff_table_hits);
			}
			std::fill(occ.begin(), occ.end(), false);
			occ[tid] = true;
			n_diff_table_hits = 1;
			len = 0;
		}
		last_pos = pos;
		// push the next match from this bucket
		VectorSeqPos& bucket = ref.hash_tables[tid].buckets_data_vectors[r->ref_bucket_id_matches_by_table[tid]];
		if(bucket_indices[tid] < bucket.size()) {
			q.push(std::make_pair(bucket[bucket_indices[tid]], tid));
			bucket_indices[tid]++;
			printf("Pushed match (%u, %u) \n", bucket[bucket_indices[tid]], tid);
		}
	}
	// add the last position
	assert(n_diff_table_hits > 0);
	assert(n_diff_table_hits-1 < params->n_tables);
	if(n_diff_table_hits >= params->min_n_hits) {
		if(n_diff_table_hits > n_best_hits) { // if more hits than best so far
			n_best_hits = n_diff_table_hits;
			if(r->ref_matches[n_diff_table_hits-1].size() < params->max_best_hits) {
				ref_match_t rm(last_pos, len);
				r->ref_matches[n_diff_table_hits-1].push_back(rm);
			}
		} else {
			// sub-optimal
			// only store if the number of hits is not too much lower than the best so far
			if(n_best_hits < params->dist_best_hit || n_diff_table_hits > (n_best_hits - params->dist_best_hit)) {
				if(r->ref_matches[n_diff_table_hits-1].size() < params->max_suboptimal_hits) {
					ref_match_t rm(last_pos, len);
					r->ref_matches[n_diff_table_hits-1].push_back(rm);
				}
			}
		}
		printf("MATCH (last) pos %u len %u n_hits %u \n", last_pos, len, n_diff_table_hits);
	}

	r->best_n_hits = (n_best_hits > 0) ? n_best_hits - 1 : 0;
	assert(r->best_n_hits < params->n_tables);
	VectorU32().swap(r->ref_bucket_id_matches_by_table); //release memory
}

void align_reads_minhash(ref_t& ref, reads_t& reads, const index_params_t* params) {
	printf("**** SRX Alignment: MinHash ****\n");

	uint32 max_windows_matched = 0;
	uint32 total_windows_matched = 0;
	uint32 diff_num_top_hits = 0;
	double start_time = omp_get_wtime();
	omp_set_num_threads(params->n_threads); // split the reads across the threads
	#pragma omp parallel reduction(+:total_windows_matched, diff_num_top_hits) reduction(max:max_windows_matched)
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
			r->ref_bucket_id_matches_by_table.resize(params->n_tables);
			for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
				minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes, params->sketch_proj_indices, t*params->sketch_proj_len);
				uint32 bucket_index = ref.hash_tables[t].bucket_indices[bucket_hash];
				r->ref_bucket_id_matches_by_table[t] = bucket_index;
			}
			collect_read_hits_contigs_inplace_merge(ref, r, params);
			//printf("Collected read %u best %u \n", i, r->best_n_hits);

			// stats
			uint32 n_contigs = 0;
			for(uint32 t = 0; t < params->n_tables; t++) {
				n_contigs += r->ref_matches[t].size();
			}
			if(n_contigs > max_windows_matched) {
				max_windows_matched = n_contigs;
			}
			total_windows_matched += n_contigs;
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
	printf("Avg number of windows matched per read %.8f \n", (float) total_windows_matched/reads.reads.size());
	printf("Avg diff of top 2 hits per read %.8f \n", (float) diff_num_top_hits/reads.reads.size());
	printf("Total number of accurate hits matching top = %d \n", acc_top);
	printf("Total number of accurate hits found = %d \n", acc_hits);
	printf("Total search time: %f sec\n", end_time - start_time);

}

// aligns the indexed reads to the indexed reference
// using simhash and permutation tables
void align_reads_lsh(ref_t& ref, reads_t& reads, const index_params_t* params) {
//	printf("**** SRX Alignment: SimHash ****\n");
//
//	// 1. sort the reads by their simhash
//	clock_t t = clock();
//	sort_reads_hash(reads);
//	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
//
//	// 2. split reads into "clusters" based on their simhash
//	t = clock();
//	VectorClusters clusters;
//	cluster_sorted_reads(reads, clusters);
//	printf("Total number of read clusters = %zu \n", clusters.size());
//	printf("Total clustering time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
//
//	// 3. for each cluster simhash, find the neighbors in the reference
//	static int perm[SIMHASH_BITLEN] = { 0 };
//	for(int i = 0; i < SIMHASH_BITLEN; i++) {
//		perm[i] = i;
//	}
//	t = clock();
//	for(uint32 p = 0; p < params->p; p++) {
//		printf("Permutation %d \n", p);
//		if(p > 0) {
//			// generate a new permutation
//			shuffle(perm);
//			permute_ref(ref, perm);
//			sort_windows_hash(ref);
//			permute_reads(clusters, perm);
//		}
//		for(uint32 i = 0; i < clusters.size(); i++) {
//			if(clusters[i].acc == 1) continue;
//
//			// binary search to find the matching ref window(s)
//			find_window_match_diffk(ref, &clusters[i], params);
//			//if(p < params->p - 1) {
//				//clusters->clusters[i].best_hamd = INT_MAX;
//			//}
//			eval_cluster_hit(&clusters[i]);
//		}
//	}
//
//	int hits = 0;
//	int acc_hits = 0;
//	int matched = 0;
//	for(uint32 i = 0; i < clusters.size(); i++) {
//		hits += clusters[i].ref_matches.size();
//		if(clusters[i].ref_matches.size() == 0) continue;
//		matched++;
//		acc_hits += eval_cluster_hit(&clusters[i]);
//		if(clusters[i].acc == 0) {
//			//printf("hash = %llx \n", clusters->clusters[i].simhash);
//			//print_read(clusters->clusters[i].reads[0]);
//			//printf("best diff = %d \n", clusters->clusters[i].best_hamd);
//			//printf("best pos = %llu \n", clusters->clusters[i].ref_matches[clusters->clusters[i].best_pos]);
//		}
//	}
//	printf("Total number of clusters matched = %d \n", matched);
//	printf("Total number of hits found = %d \n", hits);
//	printf("Total number of accurate hits found = %d \n", acc_hits);
//	printf("Total search time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	//get_stats(ref, clusters);
}

//uint32_t get_msbits32(const hash_t h, const index_params_t* params) {
//	return (h >> (SIMHASH_BITLEN - params->msbits_match));
//}

// computes the idxs of the matching simhash windows
// uses binary search
// returns -1 if no matches were found
int find_window_match_diffk(ref_t& ref, cluster_t* cluster, const index_params_t* params) {
//	const uint32_t d_q = get_msbits32(cluster->simhash, params);

//	seq_t low = 0;
//	seq_t high = ref.windows.size() - 1;
//	seq_t idx = -1;
//	while(high >= low) {
//		const seq_t mid = (low + high) / 2;
//		const uint32_t d_win = get_msbits32(ref->windows[mid].simhash, params);
//		if(d_win == d_q) {
//			idx = mid;
//			break;
//		} else if (ref->windows[mid].simhash < cluster->simhash) {
//			low = mid + 1;
//		} else {
//			high = mid - 1;
//		}
//		if(mid == 0) {
//			return -1;
//		}
//	}
//	if(idx == (seq_t) -1) {
//		return -1;
//	}
//	// find all the matches that have the same d msbits
//	seq_t l = idx;
//	if(idx > 0) {
//		l = idx - 1;
//		while (1) {
//			uint32_t d_win = get_msbits32(ref->windows[l].simhash, params);
//			if(d_win != d_q) {
//				break;
//			}
//			if(l == 0) {
//				break;
//			}
//			l--;
//		}
//		l++;
//	}
//
//	seq_t h = idx + 1;
//	while (h < ref->num_windows) {
//		uint32_t d_win = get_msbits32(ref->windows[h].simhash, params);
//		if(d_win != d_q) {
//			break;
//		}
//		h++;
//	}
//
//	//printf("Range %llu %llu \n", l, h);
//	if(cluster->ref_matches == NULL) {
//		cluster->alloc_matches = h-l+1;
//		cluster->ref_matches = (seq_t*) malloc(cluster->alloc_matches*sizeof(seq_t));
//	}
//	for(seq_t idx = l; idx < h; idx++) {
//		// check the hamming distance
//		int hammd = hamming_dist(ref->windows[idx].simhash, cluster->simhash);
//		if((hammd <= params->max_hammd)) {// && (hammd <= cluster->best_hamd)) {
//			if(cluster->num_matches == cluster->alloc_matches && cluster->alloc_matches != 0) {
//				cluster->alloc_matches <<= 1;
//				cluster->ref_matches = (seq_t*) realloc(cluster->ref_matches, cluster->alloc_matches*sizeof(seq_t));
//				if(cluster->ref_matches == NULL) {
//					printf("Could not allocate memory for the matches\n");
//					return -1;
//				}
//			}
//			if(hammd < cluster->best_hamd) {
//				cluster->best_hamd = hammd;
//				cluster->best_pos = cluster->num_matches;
//			}
//			cluster->ref_matches[cluster->num_matches] = ref->windows[idx].pos;
//			cluster->num_matches++;
//		}
//	}
	return 0;
}

void find_windows_exact(ref_t& ref, read_t* r, const index_params_t* params) {
//	seq_t low = 0;
//	seq_t high = ref->num_windows - 1;
//	seq_t idx = -1;
//	while(high >= low) {
//		seq_t mid = (low + high) / 2;
//		if(ref->windows[mid].simhash == r->simhash) {
//			idx = mid;
//			break;
//		} else if (ref->windows[mid].simhash < r->simhash) {
//			low = mid + 1;
//		} else {
//			high = mid - 1;
//		}
//		if(mid == 0) return;
//	}
//	if(idx == -1) return;
//
//	// find all the windows with the same simhash
//	seq_t l;
//	if(idx != 0) {
//		l = idx - 1;
//	} else {
//		l = 0;
//	}
//	while (idx != 0 && l >= 0) {
//		if(ref->windows[l].simhash != r->simhash) break;
//		if(l == 0) break;
//		l--;
//	}
//	if(idx != 0) l++;
//	seq_t h = idx + 1;
//	while (h < ref->num_windows) {
//		if(ref->windows[h].simhash != r->simhash) break;
//		h++;
//	}
//	h--;
//
//	if(r->ref_matches == NULL) {
//		r->alloc_matches = h-l+1;
//		r->ref_matches = (seq_t*) malloc(r->alloc_matches*sizeof(seq_t));
//	}
//	for(seq_t idx = l; idx <= h; idx++) {
//		if(r->num_matches == r->alloc_matches) {
//			r->alloc_matches <<= 1;
//			r->ref_matches = (seq_t*) realloc(r->ref_matches, r->alloc_matches*sizeof(seq_t));
//			if(r->ref_matches == NULL) {
//				printf("Could not allocate memory for the matches\n");
//				return;
//			}
//		}
//		r->ref_matches[r->num_matches] = ref->windows[idx].pos;
//		r->num_matches++;
//	}
}

char get_bucket_bits(hash_t h, int bucket) {
	char mask = 0xFF;
	return h & mask;
}

void find_windows_exact_bucket(ref_t* ref, read_t* r, const index_params_t* params, int bucket) {
//	seq_t low = 0;
//	seq_t high = ref->num_windows - 1;
//	seq_t idx = -1;
//	while(high >= low) {
//		seq_t mid = (low + high) / 2;
//		char wb = get_bucket_bits(ref->windows[mid].simhash, bucket);
//		char rb = get_bucket_bits(r->simhash, bucket);
//		if(wb == rb) {
//			idx = mid;
//			break;
//		} else if (wb < rb) {
//			low = mid + 1;
//		} else {
//			high = mid - 1;
//		}
//		if(mid == 0) return;
//	}
//	if(idx == -1) return;
//
//	// find all the windows with the same simhash
//	seq_t l;
//	if(idx != 0) {
//		l = idx - 1;
//	} else {
//		l = 0;
//	}
//	while (idx != 0 && l >= 0) {
//		if(get_bucket_bits(ref->windows[l].simhash, bucket) != get_bucket_bits(r->simhash, bucket)) break;
//		if(l == 0) break;
//		l--;
//	}
//	if(idx != 0) l++;
//	seq_t h = idx + 1;
//	while (h < ref->num_windows) {
//		if(get_bucket_bits(ref->windows[h].simhash, bucket) != get_bucket_bits(r->simhash, bucket)) break;
//		h++;
//	}
//	h--;
//
//	if(r->ref_matches == NULL) {
//		r->alloc_matches = h-l+1;
//		r->ref_matches = (seq_t*) malloc(r->alloc_matches*sizeof(seq_t));
//	}
//	for(seq_t idx = l; idx <= h; idx++) {
//		if(r->num_matches == r->alloc_matches) {
//			r->alloc_matches <<= 1;
//			r->ref_matches = (seq_t*) realloc(r->ref_matches, r->alloc_matches*sizeof(seq_t));
//			if(r->ref_matches == NULL) {
//				printf("Could not allocate memory for the matches\n");
//				return;
//			}
//		}
//		r->ref_matches[r->num_matches] = ref->windows[idx].pos;
//		r->num_matches++;
//	}
}

/* ---- Shifting/Permutation ----- */

void shift_bucket_ref(ref_t* ref, int bucket) {
//	for(seq_t i = 0; i < ref->num_windows; i++) {
//		ref->windows[i].simhash >>= (bucket*MINHASH_BUCKET_SIZE);
//	}
}

void shift_bucket_reads(reads_t* reads, int bucket) {
//	for(seq_t i = 0; i < reads->count; i++) {
//		reads->reads[i].simhash >>= (bucket*MINHASH_BUCKET_SIZE);
//	}
}

// permutes a 64-bit integer
// perm is a random permutation of integers 0..63
void perm64(uint64_t* n, int* perm) {
//	uint64_t p = 0;
//	for(int i = 0; i < SIMHASH_BITLEN; i++) {
//		int idx = perm[i];
//		p |= (((*n >> idx) & 1) << i);
//	}
//	*n = p;
}

void permute_ref(ref_t& ref, int* perm) {
//	for(seq_t i = 0; i < ref->num_windows; i++) {
//		perm64(&(ref->windows[i].simhash), perm);
//	}
}

void permute_reads(VectorClusters& clusters, int* perm) {
//	for(seq_t i = 0; i < clusters->num_clusters; i++) {
//		perm64(&(clusters->clusters[i].simhash), perm);
//	}
}

/* random integer from 0 to n-1 */
int irand(int n) {
	int r, rand_max = RAND_MAX - (RAND_MAX % n);
	/* reroll until r falls in a range that can be evenly
	 * distributed in n bins.  Unless n is comparable to
	 * to RAND_MAX, it's not *that* important really. */
	while ((r = rand()) >= rand_max);
	return r / (rand_max / n);
}

void shuffle(int *perm) {
//	int tmp;
//	int len = SIMHASH_BITLEN;
//	while(len) {
//		int j = irand(len);
//		if (j != len - 1) {
//			tmp = perm[j];
//			perm[j] = perm[len-1];
//			perm[len-1] = tmp;
//		}
//		len--;
//	}
}

/* ---- Hit Evaluation ----- */

// check how many reads in this cluster match the window positions
int eval_cluster_hit(cluster_t* cluster) {
    int matched = 0;
    for(uint32 i = 0; i < cluster->reads.size(); i++) {
        read_t r = *cluster->reads[i];
        unsigned int pos_l, pos_r;
        int strand;
        parse_read_mapping(r.name.c_str(), &pos_l, &pos_r, &strand);
        //printf("lpos %llu rpos %llu \n", pos_l, pos_r);

	int found = 0;
        for(seq_t j = pos_l - 10; j <= pos_r + 10; j++) {
        	for(seq_t idx = 0; idx < cluster->ref_matches.size(); idx++) {
            	seq_t hit_pos = cluster->ref_matches[idx];
            	if(hit_pos == j) {
            		cluster->acc = 1;
	                matched++;
	                found = 1;
	                break;
	            }
            }
        	if(found == 1) {
        		break;
        	}
        }
    }
    return matched;
}

