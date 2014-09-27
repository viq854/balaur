#ifndef INDEX_H_
#define INDEX_H_
#include "io.h"
#include "mt64.h"

#include "../rollinghashcpp/cyclichash.h"

typedef enum {SIMH, MINH, SAMPLE} algorithm;
typedef enum {OVERLAP, NON_OVERLAP, SPARSE} kmer_selection;


// universal hash function:
// single value: a*x mod n_buckets
// vector value: sum(a[i]*x[i]) mod n_buckets
struct rand_hash_function_t {
	const static uint32 w = 32; // bits per word
	minhash_t a; // odd a < 2^w
	std::vector<uint64> a_vec;
	minhash_t M; // n_butckets = 2^M

	// single int hashing
	rand_hash_function_t() {
		a = rand() + 1;
		M = w;
	}

	// vector hashing
	rand_hash_function_t(minhash_t _M, const uint32 vec_len) {
		a = 0;
		for(uint32 i = 0; i < vec_len; i++) {
			a_vec.push_back(genrand64_int64() + 1);
		}
		M = _M;
	}

	minhash_t apply(minhash_t x) const {
		return (minhash_t) a*x >> (w - M);
	}

	minhash_t apply_vector(const VectorMinHash& x, const VectorU32& indices, const uint32 vec_offset) const {
		uint64 s = 0;
		for(uint32 i = 0; i < a_vec.size(); i++) {
			s += a_vec[i]*x[indices[vec_offset + i]];
		}
		return (minhash_t) s >> (w - M);
	}
};
typedef std::vector<rand_hash_function_t> VectorHashFunctions;

typedef struct {	
	algorithm alg; 					// LSH scheme to use
	kmer_selection kmer_type; 		// scheme for extracting the kmer features

	// sequence hashing parameters => to produce sequence sketches
	uint32 ref_window_size;			// length of the reference windows to hash
	uint32 k; 						// length of the sequence kmers
	uint32 kmer_dist;				// shift between consecutive kmers
	CyclicHash* kmer_hasher;			// function used to generate kmer hashes for the sequence set
	uint32 h; 						// number of hash functions for min-hash skethes
	VectorHashFunctions minhash_functions;	// hash functions for min-hash

	// sequence kmer filtering
	uint64 max_count;				// upper bound on kmer occurrence in the reference
	uint64 min_count;				// lower bound on kmer occurrence in the read set

	// min-hash sketch hashing
	uint32 n_tables; 				// number of hash tables for the sketches
	uint32 sketch_proj_len;			// length of the sketch projection
	VectorU32 sketch_proj_indices;	// indices into the sketch for the sparse projections
	uint32 n_buckets_pow2;  		// n_buckets in a hash table = 2^n_buckets_pow2
	uint32 bucket_size;				// max number of entries to keep per bucket
	rand_hash_function_t sketch_proj_hash_func; // hash function for sketch projection vector hashing
	uint32 bucket_entry_coverage;
	uint32 contig_gap;

	// sim-hash mapping parameters
	uint32 p; 					// number of permutation tables
	uint32 msbits_match;			// number of most significant bits to match
	uint32 max_hammd; 			// maximum hamming distance to

	// multi-threading
	uint32_t n_threads;
	
	// io
	std::string in_index_fname;
	std::string out_index_fname;

	void set_default_index_params() {
		kmer_type = OVERLAP;
		ref_window_size = 150;
		k = 16;
		kmer_dist = 1;
		h = 64;
		max_count = 200;
		min_count = 0;

		n_tables = 1;
		sketch_proj_len = 4;
		n_buckets_pow2 = 16;
		bucket_size = 1000;
		bucket_entry_coverage = 10;
		contig_gap = 100;

		p = 1;
		msbits_match = 24;
		max_hammd = 10;

		n_threads = 1;
	}

} index_params_t;


typedef struct {
	std::string freq_kmer_hist_fname;
	std::string hash_func_fname;
	std::string sparse_ind_fname;

	void prep_index_files(std::string& fname) {
		freq_kmer_hist_fname += fname + std::string(".freq");
		hash_func_fname += fname + std::string(".hash");
		sparse_ind_fname += fname + std::string(".sparse");
	}
} index_files_t;


void generate_ref_windows(ref_t& ref, index_params_t* params);
void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& refidx);
void index_reads_lsh(const char* readsFname, ref_t& ref, index_params_t* params, reads_t& ridx);
void index_ref_table_i(ref_t& ref, const index_params_t* params, const seq_t i);
void index_reads_table_i(reads_t& reads, const index_params_t* params, const seq_t i);

#endif /*INDEX_H_*/
