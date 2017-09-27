#ifndef INDEX_H_
#define INDEX_H_

#include <istream>
#include <sstream>
#include <iostream>

#include "utils.h"
#include "seq.h"

#pragma once

typedef enum {SIMH, MINH, SAMPLE} algorithm;
typedef enum {OVERLAP, NON_OVERLAP, SPARSE} kmer_selection;
typedef enum {SHA1_E = 0, CITY_HASH64 = 1, PACK64 = 2} kmer_hash_alg;

#include "hash.h"

#define DISK_SYNC_PARTIAL_TABLES 0

// program parameters
typedef struct {	
	algorithm alg; 					// LSH scheme to use
	
	// LSH parameters
	kmer_selection kmer_type; 		// scheme for extracting the kmer features
	uint32 k; 						// length of the sequence kmers
	uint32 kmer_dist;			// shift between consecutive kmers
	uint32 h; 						// number of hash functions for min-hash sketches
	uint32 n_tables; 			// number of hash tables for the sketches
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

	// candidate contigs
	uint32 min_n_hits;
	uint32 dist_best_hit; 			// how many fewer than best table hits to still keep
	bool load_mhi;
	std::string precomp_contig_file_name;
	uint32 max_matched_contig_len;
	
	// voting
	uint32 k2; 							// length of the sequence kmers for vote counting
	bool precomp_k2;				// precompute k2 kmers for the reference
	uint32 delta_inlier;
	kmer_hash_alg kmer_hashing_alg;
	bool vanilla;
	bool monolith;
	int bin_size;
	int* bin_shuffle;
	int batch_size;
	bool mask_repeat_nbrs;
	int proc_contigs_thr;
	int sampling_intv;

	// alignment evaluation
	uint32 delta_x;
	uint32 n_init_anchors;
	int votes_cutoff;
	bool enable_scale;
	int mapq_scale_x;

	// multi-threading
	uint32_t n_threads;

	void set_default_index_params() {
		load_mhi = true;
		kmer_type = OVERLAP;
		h = 128;
		n_tables = 78;
		sketch_proj_len = 2;
		n_buckets_pow2 = 18;
		bucket_size = 200;
		k = 16;
		kmer_dist = 1;
		bucket_entry_coverage = 10;
		ref_window_size = 150;
		max_count = 800;
		min_count = 0;
		max_matched_contig_len = 100000;
		
		k2 = 20;
		precomp_k2 = true;
		min_n_hits = 1;
		dist_best_hit = 20;
		n_init_anchors = 10;
		delta_inlier = 10;
		delta_x = 30;
		votes_cutoff = 0;
		enable_scale = true;
		mapq_scale_x = 50;
		sampling_intv = 1;
		n_threads = 1;
		batch_size = 1;
		bin_size = 1;
		vanilla = false;
		mask_repeat_nbrs = false;
		proc_contigs_thr = 10000;//0;
		//if(ref_window_size > 150) proc_contigs_thr = 500;
		if(ref_window_size > 350) proc_contigs_thr = 20;
		kmer_hashing_alg = SHA1_E;
	}
	
	bool bin_sampling() {
		return bin_size > 1;
	}
	
	void set_bin_shuffle() {
		bin_shuffle = new int[bin_size];
		for(int i = 0; i < bin_size; i++) {
			bin_shuffle[i] = i;
		}
		shuffle(bin_shuffle, bin_size);
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
extern index_params_t* params;

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
	void release() {
		std::vector<loc_t>().swap(buckets_data);
		std::vector<uint64>().swap(bucket_offsets);
	}

} static_index_t;

// reference genome index
typedef struct {
	std::string seq; 					// reference sequence
	seq_t len;							// reference sequence length
	VectorU32 subsequence_offsets;
	std::vector<std::string> seq_names;

	MapKmerCounts kmer_hist;			// kmer occurrence histogram
	MarisaTrie high_freq_kmer_trie;		// frequent reference kmers TRIE
	VectorBool high_freq_kmer_bitmap;	// frequent reference kmers bitmap
	VectorBool ignore_kmer_bitmask;
	VectorBool ignore_window_bitmask;

	// lsh
	mutable_index_t mutable_index;
	static_index_t index;

	// voting
	std::vector<uint64> packed_32bp_kmers;
	std::vector<kmer_cipher_t> precomputed_kmer2_hashes;
	std::vector<uint16_t> precomputed_neighbor_repeats;
	std::vector<char> contig_mask;

	//std::vector<char> precomputed_local_repeats;
	//std::unordered_set<uint32> repeats;
	//std::vector<kmer_cipher_t> repeats_vec;

	void free_repeats() {
		std::vector<uint16_t>().swap(precomputed_neighbor_repeats);
	}

} ref_t;

// **** Read Set Index ****

struct ref_match_t {
	uint32_t pos;
	uint32 len;
	bool rc;
	int n_diff_bucket_hits;
	bool valid;
	ref_match_t() : pos(0), len(0), rc(0), n_diff_bucket_hits(0), valid(false) {};
	ref_match_t(seq_t _pos, uint32 _len, bool _rc, int _buckets) : pos(_pos), len(_len), rc(_rc), n_diff_bucket_hits(_buckets), valid(false) {}
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
	std::string name; 			// read name
	std::string seq;				// read sequence
	std::string rc;					// reverse complement sequence
	std::string qual;				// quality scores
	uint32 rid;

	// LSH sketches
	VectorMinHash minhashes_f;		// minhash vector
	VectorMinHash minhashes_rc;	// minhash vector for the reverse complement

	// alignment information
	std::vector<std::pair<uint64, minhash_t>> ref_bucket_matches_by_table_f;
	std::vector<std::pair<uint64, minhash_t>> ref_bucket_matches_by_table_rc;
	std::vector<ref_match_t> ref_matches;
	std::vector<bool> repeat_mask;
	kmer_cipher_t* hashes_f;
	kmer_cipher_t* hashes_rc;
	
	aln_t top_aln;
	aln_t second_best_aln;
	
	//char ref_strand;
	int n_match_f;
	char valid_minhash_f;
	char valid_minhash_rc;
	int best_n_bucket_hits;
	int true_n_bucket_hits;
	bool any_bucket_hits;
	int max_total_votes;
	int max_total_votes_low_anchors;
	char acc; 
	char top_hit_acc;
	char dp_hit_acc;
	bool collected_true_hit;
	bool processed_true_hit;
	int bucketed_true_hit;
	int comp_votes_hit;
	int n_proc_contigs;
	int strand;
	unsigned int seq_id;
	uint32_t ref_pos_l;
	uint32_t ref_pos_r;

	void free_minhash_state() {
		std::vector<minhash_t>().swap(minhashes_f);
		std::vector<minhash_t>().swap(minhashes_rc);	
	}

	read_t():  n_match_f(0),
	 valid_minhash_f(0),
	 valid_minhash_rc(0),
	 best_n_bucket_hits(0),
	 true_n_bucket_hits(0),
	 any_bucket_hits(0),
	 max_total_votes(0),
	 max_total_votes_low_anchors(0),
	 acc(0),
	 top_hit_acc(0),
	 dp_hit_acc(0),
	 collected_true_hit(0),
	 processed_true_hit(0),
	 bucketed_true_hit(0),
	 comp_votes_hit(0),
	 n_proc_contigs(0),
	 strand(0),
	seq_id(0),
	ref_pos_l(0),
	ref_pos_r(0)
	{
		top_aln.inlier_votes = 0;
        	second_best_aln.inlier_votes = 0;
		top_aln.ref_start = 0;
		second_best_aln.ref_start = 0;
		top_aln.score = 0;
		second_best_aln.score = 0;
		hashes_f = 0;
		hashes_rc = 0;
	}

	void free() {
		delete(hashes_f);
        	delete(hashes_rc);
	}
	
	// assumes that reads were generated with wgsim
	void parse_read_mapping() {
		std::istringstream is((std::string(name)));
		std::string _seqid;
		std::string refl;    
		std::string refr;
		std::getline(is, _seqid, '_');
		std::getline(is, refl, '_');
		std::getline(is, refr, '_');
		if(_seqid.compare("X") == 0) {
			seq_id = 23;
		} else {
			seq_id = atoi(_seqid.c_str());
		}
		ref_pos_l = atoi(refl.c_str());
		ref_pos_r = atoi(refr.c_str());
	}

	void get_sim_read_info(const ref_t& ref) {
#if(SIM_EVAL)
		parse_read_mapping();
		seq_id = seq_id - 1;
		if(ref.subsequence_offsets.size() > 1) {
			ref_pos_l += ref.subsequence_offsets[seq_id]; // convert to global id
			ref_pos_r += ref.subsequence_offsets[seq_id];
		}
#endif
	}

	/*void print_read() {
		printf("%s \n", name.c_str());
		for(uint32 i = 0; i < len; i++) {
			printf("%c", iupacChar[(int) seq[i]]);
		}
		printf("\n");
		for(uint32 i = 0; i < read->len; i++) {
			printf("%c", iupacChar[(int) read->rc[i]]);
		}
		printf("\n");
	}*/
	
	inline bool is_valid() {
		return (valid_minhash_f || valid_minhash_rc);
	}
	void compare_and_update_best_aln(int* n_votes, seq_t* pos, bool rc) {
		for(int i = 0; i < 2; i++) {
			if(n_votes[i] > top_aln.inlier_votes) {
				if(!pos_in_range(pos[i], top_aln.ref_start, params->delta_x)) {
					second_best_aln.inlier_votes = top_aln.inlier_votes;
					second_best_aln.total_votes = top_aln.total_votes;
					second_best_aln.ref_start = top_aln.ref_start;
				}
				// update best alignment
				top_aln.inlier_votes = n_votes[i];
				top_aln.ref_start = pos[i];
				top_aln.rc = rc;
			} else if(n_votes[i] > second_best_aln.inlier_votes) {
				if(!pos_in_range(pos[i], top_aln.ref_start, params->delta_x)) {
					second_best_aln.inlier_votes = n_votes[i];
					second_best_aln.ref_start = pos[i];
				}
			}
		}
	}
	
	void set_repeat_mask(const int k, const int n_nbrs) {
		if(repeat_mask.size() != 0) return;
		find_repeats(seq, k, n_nbrs, repeat_mask);
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

	void free() {
		for(uint32 i = 0; i < reads.size(); i++) {
			read_t* r = &reads[i];
			r->free();
		}
	}
	
	void free_minhash_data() {
		for(uint32 i = 0; i < reads.size(); i++) {
			read_t* r = &reads[i];
			r->free_minhash_state();
		}
	}

} reads_t;

void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& refidx);
void load_index_ref_lsh(const char* fastaFname, const index_params_t* params, ref_t& ref);
void store_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref);
void ref_kmer_fingerprint_stats(const char* fastaFname, index_params_t* params, ref_t& ref);

#endif /*INDEX_H_*/
