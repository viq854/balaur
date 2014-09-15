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

// generate all valid/informative reference windows
void generate_ref_windows(ref_t* ref, index_params_t* params) {
	ref->num_windows = 0;
	seq_t max_num_windows = ref->len - params->ref_window_size + 1;
	ref->windows = (ref_win_t*) malloc(max_num_windows*sizeof(ref_win_t));
	seq_t i;
	for(i = 0; i < max_num_windows; i++) {
		if(is_inform(&ref->seq[i], params->ref_window_size)) {
			ref_win_t* window = &ref->windows[ref->num_windows];
			window->pos = i;
			ref->num_windows++;
		} else {
			// skip until at least 1 potential valid kmer
			//i += params->k;
		}
	}
	ref->windows = (ref_win_t*) realloc(ref->windows, ref->num_windows*sizeof(ref_win_t));
}

// --- Simhash / Minhash ---

// index reference:
// - load fasta
// - generate valid windows (sliding window)
// - compute the hash of each window
// - sort
void index_ref_lsh(char* fastaFname, index_params_t* params, ref_t** ref_idx) {
	printf("**** SRX Reference Indexing ****\n");
		
	// 1. load the reference
	clock_t t = clock();
	ref_t* ref = fasta2ref(fastaFname);
	printf("Total ref loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 2. compute the frequency of each kmer and filter out windows to be discarded
	t = clock();
	if(params->kmer_type != SPARSE) {
		generate_ref_kmer_hist(ref, params);
		printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	
	// 3. compute valid reference windows
	t = clock();
	generate_ref_windows(ref, params);
	printf("Total number of valid windows: %llu\n", ref->num_windows);
	
	// 4. hash each window
	t = clock();
	params->max_count = (uint64_t) ceil(params->max_freq*ref->num_windows);
	for(seq_t i = 0; i < ref->num_windows; i++) {
		if(params->alg == SIMH) {
			switch (params->kmer_type) {
				case SPARSE: simhash_ref_sparse(ref, &ref->windows[i], params); break;
				case OVERLAP: simhash_ref_ovp(ref, &ref->windows[i], params); break;
				case NON_OVERLAP: simhash_ref_novp(ref, &ref->windows[i], params); break;
			}
		} else if (params->alg == MINH) {
			minhash_ref(ref, &ref->windows[i], params);
		}
	}
	printf("Total hashing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 6. sort
	t = clock();
	sort_windows_hash(ref);
	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	printf("Total number of distinct window hashes = %llu \n", get_num_distinct(ref));
	
	*ref_idx = ref;
}

void index_reads_lsh(char* readsFname, ref_t* ref, index_params_t* params, reads_t** reads_idx) {
	printf("**** SRX Read Indexing ****\n");

	// 1. load the reads (TODO: batch mode)
	clock_t t = clock();
	reads_t* reads = fastq2reads(readsFname);
	printf("Total read loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 2. compute the frequency of each kmer
	t = clock();
	if(params->kmer_type != SPARSE) {
		generate_reads_kmer_hist(reads, params);
		printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	// 3. compute the fingerprints of each read
	t = clock();
	params->min_count = (int) (params->min_freq*reads->count);
	for(int i = 0; i < reads->count; i++) {
		if(params->alg == SIMH) {
			switch (params->kmer_type) {
				case SPARSE: simhash_read_sparse(&reads->reads[i], reads->hist, ref->hist, params); break;
				case OVERLAP: simhash_read_ovp(&reads->reads[i], reads->hist, ref->hist, params); break;
				case NON_OVERLAP: simhash_read_novp(&reads->reads[i], reads->hist, ref->hist, params); break;
			}
		} else if (params->alg == MINH) {
			minhash_read(&reads->reads[i], reads->hist, ref->hist, params);
		}
	}
	printf("Total hashing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 4. sort the reads by their hash
	/*sort_reads_hash(reads);

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

	*reads_idx = reads;

	// free memory
	free(ref->hist);
	free(reads->hist);
}

// --- Sampling ---

// create hash table i 
// - hash each window
// - sort
void index_ref_table_i(ref_t* ref, index_params_t* params, int i) {
	// hash each window using sampling ids i
	clock_t t = clock();
	for(seq_t j = 0; j < ref->num_windows; j++) {
		sampling_hash_ref(ref, &ref->windows[j], params, i);
	}
	printf("Total hash table %d computation time: %.2f sec\n", i, (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// sort the hashes
	t = clock();
	sort_windows_hash(ref);
	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	printf("Total number of distinct window hashes = %llu \n", get_num_distinct(ref));
}

void index_reads_table_i(reads_t* reads, index_params_t* params, int i) {
	// hash each read using sampling ids i
	clock_t t = clock();
	for(int j = 0; j < reads->count; j++) {
		sampling_hash_read(&reads->reads[j], params, i);
	}
	printf("Total read hash table %d computation time: %.2f sec\n", i, (float)(clock() - t) / CLOCKS_PER_SEC);
}
