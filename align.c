#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "index.h"
#include "io.h"
#include "hash.h"
#include "cluster.h"

int find_window_match_diffk(ref_t* ref, cluster_t* cluster, index_params_t* params);
int find_window_match(ref_t* ref, cluster_t* cluster, index_params_t* params);
int eval_hit(cluster_t* cluster);

// aligns the indexed reads to the iindexed reference
void align_reads(ref_t* ref, reads_t* reads, index_params_t* params) {
	printf("**** SRX Alignment ****\n");
	
	// 1. sort the ref windows and the reads by their simhash
	clock_t t = clock();
	sort_windows_simhash(ref);
	sort_reads_simhash(reads);
	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	printf("Total number of distinct window hashes = %llu \n", get_num_distinct(ref));
	
	// 2. split reads into "clusters" based on their simhash 
	t = clock();
	clusters_t* clusters;
	cluster_sorted_reads(reads, &clusters);
	printf("Total number of read clusters = %d \n", clusters->num_clusters);
	printf("Total clustering time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 3. for each cluster simhash, find the neighbors in the reference
	t = clock();
	int hits = 0;
	int acc_hits = 0;
	for(int i = 0; i < clusters->num_clusters; i++) {
        // binary search to find the matching ref window(s) 
        int r = find_window_match(ref, &clusters->clusters[i], params);
        if(r < 0) {
            continue; // no match found
        }
        hits++;
        acc_hits += eval_hit(&clusters->clusters[i]);
    }
	printf("Total number of hits found = %d \n", hits);
	printf("Total number of accurate hits found = %d \n", acc_hits);
	printf("Total search time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}

// returns the idxs of the matching simhash window
// uses binary search
#define D_MSBITS_MATCH	24

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
	}
	if(idx == -1) {
		return -1;
	}
	// find all the matches that have the same d msbits
	seq_t l = idx-1;
	while (l >= 0) {
		uint32_t d_win = get_msbits32(ref->windows[l].simhash);
		uint32_t d_q = get_msbits32(cluster->simhash);
		if(d_win != d_q) {
			break;
		}
		l--;
	}
	l++;
	seq_t h = idx+1;
	while (h < ref->num_windows) {
		uint32_t d_win = get_msbits32(ref->windows[h].simhash);
		uint32_t d_q = get_msbits32(cluster->simhash);
		if(d_win != d_q) {
			break;
		}
		h++;
	}
	h--;
	cluster->ref_matches = (seq_t*) malloc((h-l+1)*sizeof(seq_t));
	for(seq_t idx = l; idx <= h; idx++) {
		// check the hamming distance
		if(hamming_dist(ref->windows[idx].simhash, cluster->simhash) <= params->max_hammd) {
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

        for(seq_t j = pos_l - 5; j <= pos_r + 5; j++) {
            int found = 0;
        	for(seq_t idx = 0; idx < cluster->num_matches; idx++) {
            	seq_t hit_pos = cluster->ref_matches[idx];
            	if(hit_pos == j) {
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