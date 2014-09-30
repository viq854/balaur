#ifndef INDEX_H_
#define INDEX_H_

#pragma once

typedef enum {SIMH, MINH, SAMPLE} algorithm;
typedef enum {OVERLAP, NON_OVERLAP, SPARSE} kmer_selection;

#include <omp.h>
#include "hash.h"

// program parameters
typedef struct {	
	algorithm alg; 					// LSH scheme to use
	kmer_selection kmer_type; 		// scheme for extracting the kmer features

	// sequence hashing parameters => to produce sequence sketches
	uint32 ref_window_size;			// length of the reference windows to hash
	uint32 k; 						// length of the sequence kmers
	uint32 kmer_dist;				// shift between consecutive kmers
	CyclicHash* kmer_hasher;			// function used to generate kmer hashes for the sequence set
	uint32 h; 						// number of hash functions for min-hash sketches
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

	// alignment evaluation
	uint32 min_n_hits;
	uint32 dist_best_hit; // how many fewer than best table hits to still keep
	uint32 max_best_hits;
	uint32 max_suboptimal_hits;
	uint32 hit_collection_interval;

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

		min_n_hits = 2;
		dist_best_hit = 2;
		max_best_hits = 1000;
		max_suboptimal_hits = 500;

		hit_collection_interval = 5000000;
	}

	// set the initial kmer hash function (rolling hash)
	void set_kmer_hash_function() {
		kmer_hasher = new CyclicHash(k, 32);
	}

	// generate random vector hash function sketch buckets
	void set_minhash_sketch_hash_function() {
		sketch_proj_hash_func = rand_hash_function_t(n_buckets_pow2, sketch_proj_len);
	}

	// generate random hash functions for min-hash sketches
	void set_minhash_hash_function() {
		for(uint32 f = 0; f < h; f++) {
			minhash_functions.push_back(rand_hash_function_t());
		}
	}

	void generate_sparse_sketch_projections() {
		rand_range_generator_t rgen;
		VectorU32 idx(h); // sketch length
		sketch_proj_indices.resize(sketch_proj_len * n_tables);
		for(uint32 i = 0; i < n_tables; i++) {
			for(uint32 k = 0; k < h; k++) {
				idx[k] = k;
			}
			// pick random indices from the sketch
			const int32_t offset = i*sketch_proj_len;
			const int32_t start = rgen.rand_in_range(h);
			sketch_proj_indices[offset] = start;
			uint32 cnt = 0;
			uint32 len = h;
			while(cnt < sketch_proj_len) {
				int j = rgen.rand_in_range(len); // exclude 0
				sketch_proj_indices[offset + cnt] = idx[j];
				idx[j] = idx[len-1];
				cnt++;
				len--;
			}
		}
	}

} index_params_t;

// **** Reference Index ****
typedef std::vector<omp_lock_t> VectorLocks;
typedef std::vector<VectorSeqPos> VectorBuckets;
typedef std::vector<VectorU32> VectorBucketIndices;

// min-hash signature index
struct buckets_t {
	uint32 n_buckets;
	uint32 next_free_bucket_index;
	VectorU32 bucket_indices;
	VectorLocks bucket_index_locks;

	VectorBuckets buckets_data_vectors;
	VectorU32 bucket_data_consumed_indices;
	omp_lock_t lock;
};
typedef std::vector<buckets_t> VectorBucketTables;

// reference genome index
typedef struct {
	std::string seq; 					// reference sequence data
	seq_t len;							// reference sequence length

	MapKmerCounts kmer_hist;			// kmer occurrence histogram
	marisa::Trie high_freq_kmer_trie;
	VectorBool ignore_kmer_bitmask;
	VectorBool ignore_window_bitmask;

	VectorBucketTables hash_tables;		// LSH min-hash index
} ref_t;

// **** Read Set Index ****

struct ref_match_t {
	const seq_t pos;
	const uint32 len;
	ref_match_t(seq_t _pos, uint32 _len) : pos(_pos), len(_len) {}
};
typedef std::vector<ref_match_t> VectorRefMatches;

typedef struct {
	uint32_t len; 					// read length
	std::string name; 				// read name
	std::string seq;				// read sequence
	std::string rc;					// reverse complement sequence
	std::string qual;				// quality scores

	// LSH fingerprints
	VectorMinHash minhashes;		// minhash vector
	char valid_minhash;

	// original mapping information
	int strand;
	uint32_t ref_pos_l;
	uint32_t ref_pos_r;

	// found ref match positions
	VectorU32 ref_bucket_id_matches_by_table;
	std::vector<VectorRefMatches> ref_matches;
	uint32 best_n_hits;

	char acc; // DEBUG: whether read matched accurately
	char top_hit_acc;

} read_t;
typedef std::vector<read_t> VectorReads;
typedef std::vector<read_t*> VectorPReads;


// collection of reads
typedef struct {
	VectorReads reads;				// read data
	MapKmerCounts kmer_hist;		// kmer histogram
	MapKmerCounts low_freq_kmer_hist;
} reads_t;


//typedef struct {
//	std::string freq_kmer_hist_fname;
//	std::string hash_func_fname;
//	std::string sparse_ind_fname;
//
//	void prep_index_files(std::string& fname) {
//		freq_kmer_hist_fname += fname + std::string(".freq");
//		hash_func_fname += fname + std::string(".hash");
//		sparse_ind_fname += fname + std::string(".sparse");
//	}
//} index_files_t;

void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& refidx);
void load_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref);
void store_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref);
void index_reads_lsh(const char* readsFname, ref_t& ref, index_params_t* params, reads_t& ridx);

#endif /*INDEX_H_*/
