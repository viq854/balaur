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

// aligns the indexed reads to the iindexed reference
void align_reads(ref_t* ref, reads_t* reads) {
	printf("**** SRX Alignment ****\n");
	
	// 1. sort the ref windows and the reads by their simhash
	clock_t t = clock();
	sort_windows_simhash(ref);
	sort_reads_simhash(reads);
	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 2. split reads into "clusters" based on their simhash 
	t = clock();
	clusters_t* clusters;
	cluster_sorted_reads(reads, &clusters);
	printf("Total number of clusters = %d \n", clusters->num_clusters);
	printf("Total clustering time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 3. for each cluster simhash, find the neighbors in the reference
	t = clock();
	int hits = 0;
	for(int i = 0; i < clusters->num_clusters; i++) {
		// binary search to find the matching ref window(s) 
		seq_t idx = find_window_match(ref, clusters->clusters[i].simhash);
		if(idx < 0) {
			continue; // no match found
		}
		hits++;
	}
	printf("Total number of hits fund = %d \n", hits);
	printf("Total search time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
}