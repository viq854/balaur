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
		if(is_inform_ref_window(&ref->seq[i], params->ref_window_size)) {
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
		ref->hist_size = params->hist_size;
		ref->hist = (uint32_t*) calloc(params->hist_size, sizeof(uint32_t));
		compute_kmer_counts(ref->seq, ref->len, params, ref->hist);
		printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	
	// 3. compute valid reference windows
	t = clock();
	generate_ref_windows(ref, params);
	printf("Total number of valid windows: %llu\n", ref->num_windows);
	
	// 4. hash each window
	t = clock();
	for(seq_t i = 0; i < ref->num_windows; i++) {
		if(params->alg == SIMH) {
			ref->windows[i].simhash = simhash(ref->seq, ref->windows[i].pos, params->ref_window_size,
					NULL, ref->hist, params, 1);
		} else if (params->alg == MINH) {
			ref->windows[i].simhash = minhash(ref->seq, ref->windows[i].pos, params->ref_window_size,
					NULL, ref->hist, params, 1);
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
		reads->hist = (uint32_t*) calloc(params->hist_size, sizeof(uint32_t));
		for(uint32_t i = 0; i < reads->count; i++) { // add the contribution of each read
			read_t* r = &reads->reads[i];
			compute_kmer_counts(r->seq, r->len, params, reads->hist);
		}
		printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	// 3. compute the fingerprints of each read
	t = clock();
	for(uint32_t i = 0; i < reads->count; i++) {
		read_t* r = &reads->reads[i];
		if(params->alg == SIMH) {
			r->simhash = simhash(r->seq, 0, r->len, reads->hist, ref->hist, params, 0);
		} else if (params->alg == MINH) {
			r->simhash = minhash(r->seq, 0, r->len, reads->hist, ref->hist, params, 0);
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
void index_ref_table_i(ref_t* ref, const index_params_t* params, const seq_t i) {
	// hash each window using sampling ids i
	clock_t t = clock();
	for(seq_t j = 0; j < ref->num_windows; j++) {
		ref->windows[j].simhash = sampling_hash(ref->seq, ref->windows[j].pos, i, params);
	}
	printf("Total hash table %llu computation time: %.2f sec\n", (uint64_t) i, (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// sort the hashes
	t = clock();
	sort_windows_hash(ref);
	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	printf("Total number of distinct window hashes = %llu \n", (uint64_t) get_num_distinct(ref));
}

void index_reads_table_i(reads_t* reads, const index_params_t* params, const seq_t i) {
	// hash each read using sampling ids i
	clock_t t = clock();
	for(uint32_t j = 0; j < reads->count; j++) {
		reads->reads[j].simhash = sampling_hash(reads->reads[j].seq, 0, i, params);
	}
	printf("Total read hash table %llu computation time: %.2f sec\n", (uint64_t) i, (float)(clock() - t) / CLOCKS_PER_SEC);
}
