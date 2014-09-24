#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <limits.h>
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

int eval_read_hit(ref_t& ref, read_t* r) {
   unsigned int pos_l, pos_r;
   int strand;
   parse_read_mapping(r->name.c_str(), &pos_l, &pos_r, &strand);

   for(uint32 i = 0; i < r->ref_matches.size(); i++) {
	   for(uint32 j = 0; j < r->ref_matches[i].size(); j++) {
		   ref_match_t match = r->ref_matches[i][j];
		   if(pos_l >= match.pos - match.len - 20 && pos_l <= match.pos + 20) {
			   r->acc = 1;
			   break;
		   }
	   }
   }
   return (r->acc == 1);
}

#define GAP_LENGTH 100

void collect_read_hits(ref_t& ref, read_t* r, const index_params_t* params) {

	// collect all the hits and their table ids
	std::vector<std::pair<seq_t, uint32> > pos_tid;
	for(uint32 t = 0; t < params->n_tables; t++) { // for each table
		buckets_t* buckets = &ref.hash_tables[t];
		uint32 bucket_index = r->ref_bucket_id_matches_by_table[t];
		if(bucket_index == buckets->n_buckets) {
			continue;
		}
		VectorSeqPos& bucket = buckets->buckets_data_vectors[bucket_index];
		for(uint32 match = 0; match < bucket.size(); match++) {
			pos_tid.push_back(std::make_pair(bucket[match], t)); // TODO: limit the number of hits collected
		}
	}
	r->ref_bucket_id_matches_by_table = VectorU32(); //release memory

	// need to find hot-spot contigs and count how many times a position occurs
	// sort and reduce
	std::sort(pos_tid.begin(), pos_tid.end());

	// construct a priority heap of matches
	// ordered by number projections matched
	r->ref_matches.resize(params->n_tables);

	// due to sampling consider different near-by positions as the frequency of the locus
	// except when they were in the same bucket! - then it should count as 1
	// near-by positions will be consecutive in the sorted list
	// # diff buckets can be found by checking bucket ids

	seq_t last_pos = pos_tid[0].first;
	uint32 last_tid = pos_tid[0].second;
	uint32 len = 0;
	VectorU32 occ(params->n_tables);
	occ[last_tid] = 1;
	for(uint32 i = 1; i < pos_tid.size(); i++) {
		seq_t pos = pos_tid[i].first;
		uint32 tid = pos_tid[i].second;
		if(pos <= last_pos + GAP_LENGTH) { // look for contigs not separated by more than GAP_LEN
			if(tid != last_tid) {
				occ[tid] = 1; // mark
			}
			len += pos - last_pos;
		} else {
			int n_occ = std::count(occ.begin(), occ.end(), 1);
			r->ref_matches[n_occ].push_back(ref_match_t(last_pos, len));
			len = 0;
			std::fill(occ.begin(), occ.end(), 0);
		}
		last_pos = pos;
		last_tid = tid;
	}
}

void align_reads_minhash(ref_t& ref, reads_t& reads, const index_params_t* params) {
	printf("**** SRX Alignment: MinHash ****\n");

	uint32 max_windows_matched = 0;
	uint32 total_windows_matched = 0;
	clock_t t = clock();
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		uint32 n_matches = 0;
		//r->ref_bucket_id_matches_by_table.resize(params->n_tables);
		for(uint32 t = 0; t < params->n_tables; t++) { // search each hash table
			VectorMinHash sketch_proj(params->sketch_proj_len);
			for(uint32 p = 0; p < params->sketch_proj_len; p++) {
				sketch_proj[p] = r->minhashes[params->sketch_proj_indices[t*params->sketch_proj_len + p]];
			}
			minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(sketch_proj);

			buckets_t* buckets = &ref.hash_tables[t];
			uint32 bucket_index = buckets->bucket_indices[bucket_hash];
			r->ref_bucket_id_matches_by_table.push_back(bucket_index);
			if(bucket_index != buckets->n_buckets) {
				n_matches += buckets->buckets_data_vectors[bucket_index].size();
			}
		}

		if(n_matches > max_windows_matched) {
			max_windows_matched = n_matches;
		}
		total_windows_matched += n_matches;
	}

	int acc_hits = 0;
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		collect_read_hits(ref, &reads.reads[i], params);
		eval_read_hit(ref, &reads.reads[i]);
		acc_hits += reads.reads[i].acc;
	}

	printf("Max number of windows matched by read %u \n", max_windows_matched);
	printf("Avg number of windows matched per read %.8f \n", (float) total_windows_matched/reads.reads.size());
	printf("Total number of accurate hits found = %d \n", acc_hits);
	printf("Total search time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

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

uint32_t get_msbits32(const hash_t h, const index_params_t* params) {
	return (h >> (SIMHASH_BITLEN - params->msbits_match));
}

// computes the idxs of the matching simhash windows
// uses binary search
// returns -1 if no matches were found
int find_window_match_diffk(ref_t& ref, cluster_t* cluster, const index_params_t* params) {
	const uint32_t d_q = get_msbits32(cluster->simhash, params);

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
	uint64_t p = 0;
	for(int i = 0; i < SIMHASH_BITLEN; i++) {
		int idx = perm[i];
		p |= (((*n >> idx) & 1) << i);
	}
	*n = p;
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
	int tmp;
	int len = SIMHASH_BITLEN; 
	while(len) {		
		int j = irand(len);				
		if (j != len - 1) {			
			tmp = perm[j];
			perm[j] = perm[len-1];
			perm[len-1] = tmp;
		} 
		len--;			
	}		
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

