#ifndef INDEX_H_
#define INDEX_H_

#pragma once

typedef enum {SIMH, MINH, SAMPLE} algorithm;
typedef enum {OVERLAP, NON_OVERLAP, SPARSE} kmer_selection;

#include "hash.h"

#define DISK_SYNC_PARTIAL_TABLES 0

// program parameters
typedef struct {	
	algorithm alg; 					// LSH scheme to use

	// LSH parameters
	kmer_selection kmer_type; 		// scheme for extracting the kmer features
	uint32 k; 						// length of the sequence kmers
	uint32 kmer_dist;				// shift between consecutive kmers
	uint32 h; 						// number of hash functions for min-hash sketches
	uint32 n_tables; 				// number of hash tables for the sketches
	uint32 n_buckets;
	uint32 sketch_proj_len;			// length of the sketch projection
	VectorU32 sketch_proj_indices;	// indices into the sketch for the sparse projections
	uint32 n_buckets_pow2;  		// n_buckets in a hash table = 2^n_buckets_pow2
	uint32 bucket_size;				// max number of entries to keep per bucket
	rand_hash_function_t sketch_proj_hash_func; // hash function for sketch projection vector hashing
	VectorHashFunctions minhash_functions;	// hash functions for min-hash
	kmer_hasher_t* kmer_hasher;		// function used to generate kmer hashes for the sequence set
	uint32 ref_window_size;			// length of the reference windows to hash
	uint32 bucket_entry_coverage;

	// sequence kmer filtering
	uint64 max_count;				// upper bound on kmer occurrence in the reference
	uint64 min_count;				// lower bound on kmer occurrence in the read set

	// alignment evaluation
	uint32 k2; 						// length of the sequence kmers for vote counting
	bool precomp_k2;				// precompute k2 kmers for the reference
	uint32 min_n_hits;
	uint32 dist_best_hit; 			// how many fewer than best table hits to still keep
	uint32 max_matched_contig_len;
	uint32 delta_inlier;
	uint32 delta_x;
	uint32 n_init_anchors;
	int votes_cutoff;
	bool enable_scale;
	int mapq_scale_x;

	// simhash mapping parameters
	uint32 p; 					// number of permutation tables
	uint32 msbits_match;		// number of most significant bits to match
	uint32 max_hammd; 			// maximum hamming distance to

	// multi-threading
	uint32_t n_threads;
	
	// io
	std::string in_index_fname;
	std::string out_index_fname;

	void set_default_index_params() {
		kmer_type = OVERLAP;
		h = 64;
		n_tables = 32;
		sketch_proj_len = 2;
		n_buckets_pow2 = 16;
		bucket_size = 200;
		k = 16;
		kmer_dist = 1;
		bucket_entry_coverage = 10;
		ref_window_size = 1000;
		max_count = 300;
		min_count = 0;
		k2 = 32;
		precomp_k2 = true;
		min_n_hits = 2;
		dist_best_hit = 25;
		max_matched_contig_len = 100000;
		n_init_anchors = 10;
		delta_inlier = 10;
		delta_x = 3;
		votes_cutoff = 50;
		enable_scale = false;
		mapq_scale_x = 100;

		n_threads = 1;
	}

	// set the initial kmer hash function (rolling hash)
	void set_kmer_hash_function() {
		kmer_hasher = new kmer_hasher_t();
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
typedef std::map<uint32, seq_t> MapKmerCounts;

// min-hash signature index
struct buckets_t {
	uint32 n_entries; // total number of entries across all buckets
	std::vector<VectorSeqPos> buckets_data_vectors;

	// per thread local buckets
	std::vector<std::vector<VectorSeqPos> > per_thread_buckets_data_vectors;
	std::vector<VectorU32> per_thread_bucket_sizes;
};

typedef struct {
	std::vector<buckets_t> per_table_buckets; // buckets for T tables
} mutable_index_t;

typedef struct {
	// stores the bucket entries across all the tables
	std::vector<loc_t> buckets_data;
	// stores offsets for each bucket id
	std::vector<uint64> bucket_offsets;
} static_index_t;

// reference genome index
typedef struct {
	std::string seq; 					// reference sequence
	seq_t len;							// reference sequence length
	VectorU32 subsequence_offsets;

	MapKmerCounts kmer_hist;			// kmer occurrence histogram
	MarisaTrie high_freq_kmer_trie;		// frequent reference kmers TRIE
	VectorBool high_freq_kmer_bitmap;	// frequent reference kmers bitmap
	VectorBool ignore_kmer_bitmask;
	VectorBool ignore_window_bitmask;
	std::vector<minhash_t> precomputed_kmer2_hashes;

	mutable_index_t mutable_index;
	static_index_t index;
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
	seq_t ref_start; 	// start position in the reference sequence
	bool rc;			// reverse complement match
	int score;
	int inlier_votes;       // number of kmers supporting the aln pos
	int total_votes;		// total number of kmers that matched
};

struct read_t {
	uint32_t len; 					// read length
	std::string name; 				// read name
	std::string seq;				// read sequence
	std::string rc;					// reverse complement sequence
	std::string qual;				// quality scores
	uint32 rid;

	// kmer hashes (unique, shuffled)
	std::vector<kmer_cipher_t> kmer_ciphers_phase1_f;
	std::vector<kmer_cipher_t> kmer_ciphers_phase1_rc;

	// LSH sketches
	VectorMinHash minhashes_f;		// minhash vector
	VectorMinHash minhashes_rc;		// minhash vector for the reverse complement
	char valid_minhash_f;
	char valid_minhash_rc;

	// kmer k2 hashes
	std::vector<std::pair<kmer_cipher_t, pos_cipher_t>> kmers_f;
	std::vector<std::pair<kmer_cipher_t, pos_cipher_t>> kmers_rc;

	// alignment information
	VectorU32 ref_bucket_id_matches_by_table;
	std::vector<std::pair<uint64, minhash_t>> ref_bucket_matches_by_table_f;
	std::vector<std::pair<uint64, minhash_t>> ref_bucket_matches_by_table_rc;
	std::vector<ref_match_t> ref_matches;
	int best_n_bucket_hits;
	int true_n_bucket_hits;
	bool any_bucket_hits;
	aln_t top_aln;
	aln_t second_best_aln;
	int max_total_votes;
	int max_total_votes_low_anchors;

	// simulation alignment info/stats
	char acc; // DEBUG: whether read matched accurately
	char top_hit_acc;
	char dp_hit_acc;
	bool collected_true_hit;
	bool processed_true_hit;
	int bucketed_true_hit;
	int comp_votes_hit;
	int n_proc_contigs;

	// original mapping information from simulations
	int strand;
	unsigned int seq_id;
	uint32_t ref_pos_l;
	uint32_t ref_pos_r;

	read_t() {
		len = 0;
		rid = 0;

		valid_minhash_f = 0;
		valid_minhash_rc = 0;
		best_n_bucket_hits = 0;
		true_n_bucket_hits = 0;
		any_bucket_hits = false;

		// alignment info
		top_aln.score = 0;
		top_aln.ref_start = 0;
		top_aln.inlier_votes = 0;
		top_aln.total_votes = 0;
		second_best_aln.inlier_votes = 0;
		second_best_aln.total_votes = 0;
		max_total_votes = 0;
		max_total_votes_low_anchors = 0;

		// simulation alignment info/stats
		acc = 0;
		collected_true_hit = 0;
		processed_true_hit = false;
		bucketed_true_hit = 0;
		comp_votes_hit = 0;
		n_proc_contigs = 0;
		dp_hit_acc = 0;
		top_hit_acc = 0;
		strand = 0;
		seq_id = 0;
		ref_pos_l = 0;
		ref_pos_r = 0;
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

void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& refidx);
void load_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref);
void store_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref);
void index_reads_lsh(const char* readsFname, ref_t& ref, index_params_t* params, reads_t& ridx);
void ref_kmer_fingerprint_stats(const char* fastaFname, index_params_t* params, ref_t& ref);

#endif /*INDEX_H_*/
