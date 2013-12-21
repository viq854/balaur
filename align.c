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

// returns the idx of the matching simhash window
// uses binary search
seq_t find_window_match(ref_t* ref, simhash_t h) {
	seq_t low = 0;
	seq_t high = ref->num_windows;
	
	while(high > low) {
		seq_t mid = (low + high) / 2;
		if(ref->windows[mid].simhash == h) {
			return mid;
		} else if (ref->windows[mid].simhash < h) {
			low = mid + 1;
		} else {
			high = mid - 1;
		}
	}
	return -1;
}

// check if any read in this cluster matches the window position
int eval_hit(cluster_t* cluster, seq_t hit_pos) {
	read_t r = *cluster->reads[0];
	unsigned int pos_l, pos_r;
	int strand;
	parse_read_mapping(r.name, &pos_l, &pos_r, &strand); 

    for(int j = -10; j <= pos_r - pos_l + 1 + 10; j++) {
        if(hit_pos == (pos_l + j) - 1) {
            return 1;
        }
    }
    return 0;
}

// aligns the indexed reads to the iindexed reference
void align_reads(ref_t* ref, reads_t* reads) {
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
		seq_t idx = find_window_match(ref, clusters->clusters[i].simhash);
		if(idx < 0) {
			continue; // no match found
		}
		hits++;
		
		if(eval_hit(&clusters->clusters[i], ref->windows[idx].pos)) {
			acc_hits++;
		}
	}
	printf("Total number of hits found = %d \n", hits);
	printf("Total number of accurate hits found = %d \n", acc_hits);
	printf("Total search time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
}