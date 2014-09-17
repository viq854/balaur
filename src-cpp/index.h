#ifndef INDEX_H_
#define INDEX_H_
#include "io.h"

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
	VectorU32 sparse_kmers;			// indices into the read of the sparse kmers

	uint32_t h; 					// number of hash functions for min-hash
	VectorHash rand_hash_pads;			// hash functions for min-hash
	uint32_t n_min_matched;			// minimum number of min-hashes for a match
	
	// kmer weighing
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
	std::string in_index_fname;
	std::string out_index_fname;

	// multi-threading
	uint32_t n_threads;

	void set_default_index_params() {
		k = 16;
		kmer_dist = 1;
		m = 10;
		max_range = k + 10;
		p = 1;
		msbits_match = 24;
		h = 64;
		n_min_matched = 5;
		min_freq = 0.000001;
		max_freq = 0.6;
		ref_window_size = 100;
		max_hammd = 10;
		kmer_type = SPARSE;
		n_threads = 1;
	}

} index_params_t;


void generate_ref_windows(ref_t& ref, index_params_t* params);
void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& refidx);
void index_reads_lsh(const char* readsFname, ref_t& ref, index_params_t* params, reads_t& ridx);
void index_ref_table_i(ref_t& ref, const index_params_t* params, const seq_t i);
void index_reads_table_i(reads_t& reads, const index_params_t* params, const seq_t i);

#endif /*INDEX_H_*/
