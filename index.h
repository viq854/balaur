#ifndef INDEX_H_
#define INDEX_H_
#include "io.h"

#define DEBUG 0

typedef enum {SIMH, MINH, SAMPLE} algorithm;
typedef enum {OVERLAP, NON_OVERLAP, SPARSE} kmer_selection;

typedef struct {	
	algorithm alg; 					// LSH scheme to use
	kmer_selection kmer_type; 		// scheme for extracting the kmer features

	// hashing parameters
	uint32_t ref_window_size;		// length of the reference windows to hash
	uint32_t k; 					// length of the kmers
	uint32_t kmer_dist;				// shift between consecutive kmers (non-sparse type only)
	uint32_t m; 					// number of kmers to extract for the sparse kmers type
	uint32_t max_range; 			// range from which to locally generate sparse kmers
	uint32_t* sparse_kmers;			// indices into the read of the sparse kmers
	uint32_t h; 					// number of hash functions for min-hash
	hash_t* rand_hash_pads;			// hash functions for min-hash
	
	// kmer weighing
	uint64_t hist_size; 			// length of the kmer reference histogram
	float min_freq;
	float max_freq;
	uint64_t min_count;
	uint64_t max_count;

	// mapping parameters
	uint32_t msbits_match;			// number of most significant bits to match
	uint32_t p; 					// number of permutation tables
	uint32_t max_hammd; 			// maximum hamming distance to
	
	uint32_t s; 					// length of the hash vector

	// io
	const char* in_idx_fname;
	const char* out_idx_fname;

	// multi-threading
	uint32_t n_threads;

} index_params_t;


void generate_ref_windows(ref_t* ref, index_params_t* params);
void index_ref_lsh(char* fastaFname, index_params_t* params, ref_t** refidx);
void index_reads_lsh(char* readsFname, ref_t* ref, index_params_t* params, reads_t** ridx);
void index_ref_table_i(ref_t* ref, const index_params_t* params, const seq_t i);
void index_reads_table_i(reads_t* reads, const index_params_t* params, const seq_t i);

#endif /*INDEX_H_*/
