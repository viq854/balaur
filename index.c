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

void index_reads(char* readsFname, index_params_t* params) {
	printf("**** SRX Read Indexing ****\n");
	
	// 1. load the reads (TODO: batch mode)
	clock_t t = clock();
	reads_t* reads = fastq2reads(readsFname);
	printf("Total read loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 2. compute the fingerprints of each read
	t = clock();
	for(int i = 0; i < reads->count; i++) {
		simhash(&reads->reads[i], params);
		printf("read %d hash = %llx \n", i, reads->reads[i].simhash);
	}
	printf("Total simhash computation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 3. split reads into "clusters" based on their fp
	t = clock();
	clusters_t* clusters;
	cluster_reads(reads, &clusters);
	
	printf("Total number of clusers = %d \n", clusters->num_clusters);
	for(int i = 0; i < clusters->num_clusters; i++) {
		printf("cluster = %d, simhash = %llx, size = %d \n", i, clusters->clusters[i].simhash, clusters->clusters[i].size);
	}
	printf("Total clustering time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}