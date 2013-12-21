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
	
	// 2. compute the frequency of each kmer
	int* histogram;
	t = clock();
	generate_reads_kmer_hist(reads, params, &histogram);
	printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 3. compute the fingerprints of each read
	t = clock();
	params->min_count = (int) (params->min_freq*reads->count);
	for(int i = 0; i < reads->count; i++) {
		simhash_read(&reads->reads[i], histogram, params);
		//cityhash(&reads->reads[i]);
		//printf("read %d hash = %llx \n", i, reads->reads[i].simhash);
	}
	printf("Total simhash computation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 4. sort the reads by their simhash
	/*sort_reads_simhash(reads);
	
	// 5. split reads into "clusters" based on their simhash 
	t = clock();
	clusters_t* clusters;
	cluster_sorted_reads(reads, &clusters);
	printf("Total number of clusters = %d \n", clusters->num_clusters);
	printf("Total clustering time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	//int count = 0;
	//for(int i = 0; i < clusters->num_clusters; i++) {
		//if((count < 20) && (clusters->clusters[i].size > 1)) {
			//printf("cluster = %d, simhash = %llx, size = %d \n", i, clusters->clusters[i].simhash, clusters->clusters[i].size);
			//count++;
			//for(int j = 0; j < clusters->clusters[i].size; j++) {
				//print_read(clusters->clusters[i].reads[j]);
			//}
		//}
	//}
	
	// 6. collapse the clusters that are close to each other in Hamming distance
	t = clock();
	int num_collapsed = collapse_clusters(clusters, params);
	printf("Collapsed %d clusters \n", num_collapsed);
	printf("Total number of clusters remaining = %d \n", clusters->num_clusters - num_collapsed);
	printf("Total clustering collapse time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	*/
	
	//free_reads(reads);
	//free(histogram);
	//free(clusters->clusters); //TODO
}

void index_ref(char* fastaFname, index_params_t* params) {
	printf("**** SRX Reference Indexing ****\n");
		
	// 1. load the reference
	clock_t t = clock();
	ref_t* ref = fasta2ref(fastaFname);
	printf("Total read loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 2. compute the frequency of each kmer and filter out windows to be discarded
	int* histogram;
	t = clock();
	generate_ref_kmer_hist(ref, params, &histogram);
	printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	printf("Total number of valid windows: %llu\n", ref->num_windows);
	
	// 3. hash each window of the reference
	t = clock();
	params->max_count = (uint64_t) (params->max_freq*ref->num_windows);
	for(seq_t i = 0; i < ref->num_windows; i++) {
		simhash_ref(ref, &ref->windows[i], histogram, params);
	}
	printf("Total simhash computation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 4. sort the windows by their simhash
	//sort_windows_simhash(ref);
}