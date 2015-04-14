#ifndef INDEX_H_
#define INDEX_H_

#pragma once

typedef enum {SIMH, MINH, SAMPLE} algorithm;
typedef enum {OVERLAP, NON_OVERLAP, SPARSE} kmer_selection;

#include <omp.h>
#include "hash.h"

#define INDEX_READS_RC 1
#define INDEX_READS_REF 0
#define DISK_SYNC_PARTIAL_TABLES 0

// program parameters
typedef struct {	
	algorithm alg; 					// LSH scheme to use
	kmer_selection kmer_type; 		// scheme for extracting the kmer features

	// sequence hashing parameters => to produce sequence sketches
	uint32 ref_window_size;			// length of the reference windows to hash
	uint32 k; 						// length of the sequence kmers
	uint32 k2; 						// length of the sequence kmers for vote counting
	uint32 kmer_dist;				// shift between consecutive kmers
	CyclicHash* kmer_hasher;		// function used to generate kmer hashes for the sequence set
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

	// seed extension / chaining
	uint32 bandw;
	uint32 max_chain_gap;
	uint32 match;
	uint32 mismatch;
	uint32 gap_open;
	uint32 gap_extend;
	uint32 zdrop;
	int8_t score_matrix[25];

	// alignment evaluation
	uint32 min_n_hits;
	uint32 dist_best_hit; // how many fewer than best table hits to still keep
	uint32 max_matched_contig_len;
	uint32 n_top_buckets_search;
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
		k2 = 16;
		kmer_dist = 1;
		h = 64;
		max_count = 200;
		min_count = 0;

		n_tables = 1;
		sketch_proj_len = 4;
		n_buckets_pow2 = 16;
		bucket_size = 200;
		bucket_entry_coverage = 10;
		contig_gap = 100;

		p = 1;
		msbits_match = 24;
		max_hammd = 10;

		n_threads = 1;

		min_n_hits = 2;
		dist_best_hit = 2;
		max_matched_contig_len = 100000;
		n_top_buckets_search = 1;

		max_best_hits = 100;
		max_suboptimal_hits = 500;
		hit_collection_interval = 200000000;

		bandw = 100;
		max_chain_gap = 500;
		match = 1;
		mismatch = 4;
		gap_open = 6;
		gap_extend = 1;
		zdrop = 100;
		int i, j, k;
		for (i = k = 0; i < 4; ++i) {
			for (j = 0; j < 4; ++j) {
				score_matrix[k++] = i == j? match : -mismatch;
			}
			score_matrix[k++] = -1; // ambiguous base
		}
		for (j = 0; j < 5; ++j) {
			score_matrix[k++] = -1;
		}
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
typedef std::vector<VectorSeqPos> VectorBuckets;
typedef std::vector<std::vector<VectorSeqPos> > VectorPerThreadBuckets;
typedef std::vector<VectorU32> VectorPerThreadIndices;
typedef std::vector<VectorU32> VectorPerThreadSizes;

// min-hash signature index
struct buckets_t {
	uint32 n_buckets; // B

	// global
	uint32 next_free_bucket_index;
	VectorU32 bucket_indices;
	VectorBuckets buckets_data_vectors;

	// per thread buckets
	VectorU32 per_thread_next_free_bucket_index;
	VectorPerThreadIndices per_thread_bucket_indices;
	VectorPerThreadBuckets per_thread_buckets_data_vectors;
	VectorPerThreadSizes per_thread_bucket_sizes;
};
typedef std::vector<buckets_t> VectorBucketTables;

// reference genome index
typedef struct {
	std::string seq; 					// reference sequence
	std::string seq_RC; 				// reference RC
	seq_t len;							// reference sequence length
	VectorU32 subsequence_offsets;

	MapKmerCounts kmer_hist;			// kmer occurrence histogram
	marisa::Trie high_freq_kmer_trie;	// frequent reference kmers
	marisa::Trie high_freq_kmer_trie_RC;// frequent RC kmers
	VectorBool ignore_kmer_bitmask;
	VectorBool ignore_kmer_bitmask_RC;
	VectorBool ignore_window_bitmask;
	VectorBool ignore_window_bitmask_RC;

	VectorBucketTables hash_tables;		// LSH min-hash index, T

	std::vector<minhash_t> precomputed_kmer2_hashes;
} ref_t;

// **** Read Set Index ****

struct ref_match_t {
	uint32_t pos;
	uint32 len;
	bool rc;
	ref_match_t() : pos(0), len(0), rc(0) {};
	ref_match_t(seq_t _pos, uint32 _len, bool _rc) : pos(_pos), len(_len), rc(_rc) {}
};
typedef std::vector<ref_match_t> VectorRefMatches;

struct aln_t {
	seq_t ref_start, ref_end; 	// [rb,re): reference sequence in the alignment
	bool rc;
	int read_start, read_end;   // [qb,qe): query sequence in the alignment
	int score;
	int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
	int sub;        // 2nd best SW score
	int w;          // actual band width used in extension
};

struct read_t {
	uint32_t len; 					// read length
	std::string name; 				// read name
	std::string seq;				// read sequence
	std::string rc;					// reverse complement sequence
	std::string qual;				// quality scores

	// LSH sketches
	VectorMinHash minhashes;		// minhash vector
	VectorMinHash minhashes_rc;		// minhash vector for the reverse complement
	char valid_minhash;
	char valid_minhash_rc;

	// kmer hashes
	std::vector<std::pair<minhash_t, uint32>> kmers_f;
	std::vector<std::pair<minhash_t, uint32>> kmers_rc;

	// Mappings
	VectorU32 ref_bucket_id_matches_by_table;
	std::vector< VectorSeqPos * > ref_bucket_matches_by_table;
	std::vector<VectorRefMatches> ref_matches;
	std::vector<int> ref_match_sizes;

	std::vector<int16_t>* ref_brackets_f;
	std::vector<int16_t>* ref_brackets_rc;
	std::vector<int>* ref_brackets_dirty_f;
	std::vector<int>* ref_brackets_dirty_rc;
	uint32 rid;

	int best_n_bucket_hits;
	bool any_bucket_hits;
	int n_max_votes;
	aln_t aln;
	int max_votes;
	int max_votes_second_best;

	int max_possible_votes; //number of votes maximally possible
	int max_votes_noransac; //ignoring ransac, how many matches would there have been
	int max_votes_noransac_second_best; //ignoring ransac, how many matches would the second have had
	

	char acc; // DEBUG: whether read matched accurately
	char top_hit_acc;
	char dp_hit_acc;
	bool collected_true_hit;
	bool processed_true_hit;
	int bucketed_true_hit;
	int comp_votes_hit;

	// original mapping information from simulations
	int strand;
	unsigned int seq_id;
	uint32_t ref_pos_l;
	uint32_t ref_pos_r;

	read_t() {
		valid_minhash = 0;
		valid_minhash_rc = 0;
		best_n_bucket_hits = 0;
		any_bucket_hits = false;
		n_max_votes = 0;
		aln.score = 0;
		aln.ref_start = 0;
		max_votes = 0;
		max_votes_second_best = 0;
		acc = 0;
		collected_true_hit = 0;
		processed_true_hit = false;
		bucketed_true_hit = 0;
		comp_votes_hit = 0;
		strand = 0;
		seq_id = 0;
		ref_pos_l = 0;
		ref_pos_r = 0;
		dp_hit_acc = 0;
		top_hit_acc = 0;
		len = 0;
		rid = 0;
		ref_brackets_f = 0;
		ref_brackets_rc = 0;
		ref_brackets_dirty_f = 0;
		ref_brackets_dirty_rc = 0;

		max_possible_votes = 0;
		max_votes_noransac = 0;
		max_votes_noransac_second_best = 0;
	}
};
typedef std::vector<read_t> VectorReads;
typedef std::vector<read_t*> VectorPReads;


// collection of reads
typedef struct {
	const char* fname;
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
