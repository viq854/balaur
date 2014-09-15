#ifndef INDEX_H_
#define INDEX_H_
#include "io.h"

#define DEBUG 0

typedef enum {SIMH, MINH, SAMPLE} algorithm;
typedef enum {OVERLAP, NON_OVERLAP, SPARSE} kmer_selection;

typedef struct {	
	algorithm alg; 				// LSH scheme to use
	kmer_selection kmer_type; 	// scheme for extracting the kmer features

	// hashing parameters
	int ref_window_size;		// length of the reference windows to hash
	int k; 						// length of the kmers
	int m; 						// number of kmers to extract for the sparse kmers type
	int max_range; 				// range from which to locally generate sparse kmers
	int* sparse_kmers;			// indices into the read of the sparse kmers
	int h; 						// number of hash functions for min-hash
	hash_t* rand_hash_pads;	// hash functions for min-hash
	
	// kmer weighing
	uint64_t hist_size; 		// length of the kmer reference histogram
	float min_freq;
	float max_freq;
	uint64_t min_count;
	uint64_t max_count;

	// mapping parameters
	int p; 						// number of permutation tables
	int max_hammd; 				// maximum hamming distance to
	
	int s; // length of the hash vector

	// io
	const char* in_idx_fname;
	const char* out_idx_fname;

} index_params_t;


void generate_ref_windows(ref_t* ref, index_params_t* params);
void index_ref_lsh(char* fastaFname, index_params_t* params, ref_t** refidx);
void index_reads_lsh(char* readsFname, ref_t* ref, index_params_t* params, reads_t** ridx);
void index_ref_table_i(ref_t* ref, index_params_t* params, int i);
void index_reads_table_i(reads_t* reads, index_params_t* params, int i);

#endif /*INDEX_H_*/
