#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "index.h"
#include "io.h"
#include "hash.h"
#include "cluster.h"

int find_window_match_diffk(ref_t* ref, cluster_t* cluster, index_params_t* params);
int find_window_match(ref_t* ref, cluster_t* cluster, index_params_t* params);
int eval_hit(cluster_t* cluster);

void shuffle(int* perm);
void permute_ref(ref_t* ref, int perm[]);
void permute_reads(clusters_t* reads, int perm[]);


void diff_stats(ref_t* ref, clusters_t* clusters) {
	int diff_range = 30;
	int* diff_hist = (int*) calloc(diff_range+1, sizeof(int));
	for(int i = 0; i < clusters->num_clusters; i++) {
		cluster_t c = clusters->clusters[i];
		for(int j = 0; j < c.size; j++) {
			read_t* r = c.reads[j];
			unsigned int pos_l, pos_r;
			int strand;
			parse_read_mapping(r->name, &pos_l, &pos_r, &strand);
			unsigned int true_pos = pos_l - 1;
			for(seq_t k = 0; k < ref->num_windows; k++) {
				if(true_pos == ref->windows[k].pos) {
					int d = hamming_dist(ref->windows[k].simhash, c.simhash); 
					if(d < diff_range) {
						diff_hist[d]++;
					} else {
						diff_hist[diff_range]++;
					}
					break;
				}
			}
		}
	}
	
	for(int i = 0; i <= diff_range; i++) {
		printf("diff = %d count = %d \n", i, diff_hist[i]);
	}
	
	free(diff_hist);
}

void get_stats(ref_t* ref, clusters_t* clusters) {
	int printed = 0;
	for(int i = 0; i < clusters->num_clusters; i++) {
		cluster_t c = clusters->clusters[i];
		if(c.acc == 0) {
			if(printed > 50) break;
			printed++;
			printf("Cluster %d size = %d hash %llx \n", i, c.size, c.simhash);
			printf("best diff = %d \n", c.best_hamd);
			printf("best pos = %llu \n", c.ref_matches[c.best_pos]);
			for(int j = 0; j < c.size; j++) {
				read_t* r = c.reads[j];
				print_read(r);
				unsigned int pos_l, pos_r;
				int strand;
				parse_read_mapping(r->name, &pos_l, &pos_r, &strand);
				unsigned int true_pos = pos_l - 1;
				for(seq_t k = 0; k < ref->num_windows; k++) {
					if(true_pos == ref->windows[k].pos) {
						printf("hamm dist true = %d simhash = %llx \n", hamming_dist(ref->windows[k].simhash, c.simhash), ref->windows[k].simhash);
					}
					if(c.ref_matches[c.best_pos] == ref->windows[k].pos) {
						printf("hamm dist best = %d simhash = %llx\n", hamming_dist(ref->windows[k].simhash, c.simhash), ref->windows[k].simhash); 
					}
				}
			}
		}
	}
}

// aligns the indexed reads to the iindexed reference
void align_reads(ref_t* ref, reads_t* reads, index_params_t* params) {
	printf("**** SRX Alignment ****\n");
	
	// 1. sort the reads by their simhash
	clock_t t = clock();
	sort_reads_simhash(reads);
	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 2. split reads into "clusters" based on their simhash 
	t = clock();
	clusters_t* clusters;
	cluster_sorted_reads(reads, &clusters);
	printf("Total number of read clusters = %d \n", clusters->num_clusters);
	printf("Total clustering time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);


	//diff_stats(ref, clusters);

	// 3. for each cluster simhash, find the neighbors in the reference
	static int perm[SIMHASH_BITLEN] = { 0 };
	for(int i = 0; i < SIMHASH_BITLEN; i++) {
		perm[i] = i;
	}
	t = clock();
	for(int p = 0; p < params->p; p++) {
		printf("Permutation %d \n", p);
		if(p > 0) {
			// generate a new permutation
			shuffle(perm);
			permute_ref(ref, perm);
			sort_windows_simhash(ref);
			permute_reads(clusters, perm);
		}
		for(int i = 0; i < clusters->num_clusters; i++) {
			if(clusters->clusters[i].acc == 1) continue;
			
			// binary search to find the matching ref window(s) 
			find_window_match_diffk(ref, &clusters->clusters[i], params);
			if(p < params->p - 1) {
				//clusters->clusters[i].best_hamd = INT_MAX;
			}
			eval_hit(&clusters->clusters[i]);
		}
	}
	
	int hits = 0;
	int acc_hits = 0;
	int matched = 0;
	for(int i = 0; i < clusters->num_clusters; i++) {
		hits += clusters->clusters[i].num_matches;
		if(clusters->clusters[i].num_matches == 0) continue;
		matched++;
		acc_hits += eval_hit(&clusters->clusters[i]);
		if(clusters->clusters[i].acc == 0) {
			//printf("hash = %llx \n", clusters->clusters[i].simhash);
			//print_read(clusters->clusters[i].reads[0]);
			//printf("best diff = %d \n", clusters->clusters[i].best_hamd);
			//printf("best pos = %llu \n", clusters->clusters[i].ref_matches[clusters->clusters[i].best_pos]);
		}	
	}
	printf("Total number of clusters matched = %d \n", matched);
	printf("Total number of hits found = %d \n", hits);
	printf("Total number of accurate hits found = %d \n", acc_hits);
	printf("Total search time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	//get_stats(ref, clusters);
}

// returns the idxs of the matching simhash window
// uses binary search
#define D_MSBITS_MATCH	16

uint32_t get_msbits32(simhash_t h) {
	return (h >> (SIMHASH_BITLEN - D_MSBITS_MATCH));
}

int find_window_match_diffk(ref_t* ref, cluster_t* cluster, index_params_t* params) {
	seq_t low = 0;
	seq_t high = ref->num_windows - 1;
	seq_t idx = -1;
	while(high >= low) {
		seq_t mid = (low + high) / 2;
		uint32_t d_win = get_msbits32(ref->windows[mid].simhash);
		uint32_t d_q = get_msbits32(cluster->simhash);
		if(d_win == d_q) {
			idx = mid;
			break;
		} else if (ref->windows[mid].simhash < cluster->simhash) {
			low = mid + 1;
		} else {
			high = mid - 1;
		}
		if(mid == 0) {
			return -1;
		}
	}
	if(idx == -1) {
		return -1;
	}
	// find all the matches that have the same d msbits
	seq_t l;
	if(idx != 0) {
		l = idx - 1;
	} else {
		l = 0;
	}
	while (idx != 0 && l >= 0) {
		uint32_t d_win = get_msbits32(ref->windows[l].simhash);
		uint32_t d_q = get_msbits32(cluster->simhash);
		if(d_win != d_q) {
			break;
		}
		if(l == 0) {
			break;
		}
		l--;
	}
	if(idx != 0) l++;
	seq_t h = idx + 1;
	while (h < ref->num_windows) {
		uint32_t d_win = get_msbits32(ref->windows[h].simhash);
		uint32_t d_q = get_msbits32(cluster->simhash);
		if(d_win != d_q) {
			break;
		}
		h++;
	}
	h--;

	//printf("Range %llu %llu \n", l , h);	
	if(cluster->ref_matches == NULL) {
		cluster->alloc_matches = h-l+1;
		cluster->ref_matches = (seq_t*) malloc(cluster->alloc_matches*sizeof(seq_t));
	}
	for(seq_t idx = l; idx <= h; idx++) {
		// check the hamming distance
		int hammd = hamming_dist(ref->windows[idx].simhash, cluster->simhash);
		if((hammd <= params->max_hammd)) {// && (hammd <= cluster->best_hamd)) {
			if(cluster->num_matches == cluster->alloc_matches && cluster->alloc_matches != 0) {
				cluster->alloc_matches <<= 1;
				cluster->ref_matches = (seq_t*) realloc(cluster->ref_matches, cluster->alloc_matches*sizeof(seq_t));
				if(cluster->ref_matches == NULL) {
					printf("Could not allocate memory for the matches\n");
					return -1;
				}
			}
			if(hammd < cluster->best_hamd) {
				cluster->best_hamd = hammd;
				cluster->best_pos = cluster->num_matches;
			}
			cluster->ref_matches[cluster->num_matches] = ref->windows[idx].pos;
			cluster->num_matches++;
		}
	}
	return 0;
}

// returns the idx of the matching simhash window
// uses binary search
int find_window_match(ref_t* ref, cluster_t* cluster, index_params_t* params) {
	simhash_t h = cluster->simhash;
	seq_t low = 0;
	seq_t high = ref->num_windows - 1;
	
	while(high >= low) {
		seq_t mid = (low + high) / 2;
		if(ref->windows[mid].simhash == h) {
			cluster->ref_matches = (seq_t*) malloc(1*sizeof(seq_t));
			cluster->num_matches = 1;
			cluster->ref_matches[0] = ref->windows[mid].pos;
			return 0;
		} else if (ref->windows[mid].simhash < h) {
			low = mid + 1;
		} else {
			high = mid - 1;
		}
		if(mid == 0) {
			return -1;
		}
	}
	return -1;
}

// check how many reads in this cluster match the window positions
int eval_hit(cluster_t* cluster) {
    int matched = 0;
    for(int i = 0; i < cluster->size; i++) {
        read_t r = *cluster->reads[i];
        unsigned int pos_l, pos_r;
        int strand;
        parse_read_mapping(r.name, &pos_l, &pos_r, &strand);
        //printf("lpos %llu rpos %llu \n", pos_l, pos_r);

	int found = 0;
        for(seq_t j = pos_l - 10; j <= pos_r + 10; j++) {
        	for(seq_t idx = 0; idx < cluster->num_matches; idx++) {
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

void permute_ref(ref_t* ref, int* perm) {
	for(seq_t i = 0; i < ref->num_windows; i++) {
		perm64(&(ref->windows[i].simhash), perm);
	}
}

void permute_reads(clusters_t* clusters, int* perm) {
	for(seq_t i = 0; i < clusters->num_clusters; i++) {
		perm64(&(clusters->clusters[i].simhash), perm);
	}
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

