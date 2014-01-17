#ifndef INDEX_H_
#define INDEX_H_
#include "io.h"

typedef struct {
	int p; // number of permutation tables
	int k; // length of the k-mers
	int m; // number of k-mers for sparse k-mers
	int h; // number of hash functions for min-hash
	int* sparse_kmers;
	simhash_t* rand_hash_pads;
	
	uint64_t hist_size; // length of the k-mer freq histogram
	
	float min_freq;
	float max_freq;
	uint64_t min_count;
	uint64_t max_count;
	
	int max_hammd; // maximum hamming distance to
	
	int ref_window_size;
	
	int s; // length of the hash vector
} index_params_t;

void index_ref_simhash(char* fastaFname, index_params_t* params, ref_t** refidx);
void index_ref_windows(char* fastaFname, index_params_t* params, ref_t** refidx);
void index_ref_table_i(ref_t* ref, index_params_t* params, int i);

void index_reads_simhash(char* readsFname, ref_t* ref, index_params_t* params, reads_t** ridx);
void index_reads_table_i(reads_t* reads, index_params_t* params, int i);

#endif /*INDEX_H_*/
