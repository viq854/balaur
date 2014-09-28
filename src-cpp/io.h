#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "types.h"

#include <marisa.h>

static const unsigned char iupacChar[5] =  {'A', 'G', 'C', 'T', 'N'};
static const unsigned char nt4_complement[5] = {3/*A*/, 2/*G*/, 1/*C*/, 0/*T*/, 4/*N*/};

// encoding: A=0, G=1, C=2, T=3, N=4
static const unsigned char nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4,0/*A*/,4,2/*C*/,4,4,4,1/*G*/,4,4,4,4,4,4,4/*N*/,4,
	4, 4, 4, 4,  3/*T*/,4,4,4,4, 4, 4, 4,  4, 4, 4, 4,
	4,0/*a*/,4,2/*c*/,4,4,4,1/*g*/,4,4,4,4,4,4,4/*n*/,4,
	4, 4, 4, 4,  3/*t*/,4,4,4,4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// reference genome window
typedef struct {
	seq_t pos;
	hash_t simhash;
	VectorMinHash minhashes;
} ref_win_t;
typedef std::vector<ref_win_t*> VectorWindowPtr;
typedef std::map<seq_t, ref_win_t> MapPos2Window;
typedef std::vector<std::map<minhash_t, VectorWindowPtr> > VectorMinHashMaps;
typedef std::map<seq_t, uint32> MapPos2MinCount;

typedef std::vector<VectorSeqPos> VectorBuckets;
struct buckets_t {
	uint32 n_buckets;
	uint32 next_free_bucket_index;
	VectorU32 bucket_indices;
	//VectorU32 bucket_sizes;
	VectorBuckets buckets_data_vectors;
};
typedef std::vector<buckets_t> VectorBucketTables;
typedef std::vector<VectorU32> VectorBucketIndices;

struct minhash_matrix_t {
	std::vector<VectorMinHash> h_minhash_cols;
	uint32 oldest_col_index;
};

// reference genome
typedef struct {
	std::string seq; 					// reference sequence data
	seq_t len;							// reference sequence length
	MapKmerCounts kmer_hist;			// kmer occurrence histogram
	marisa::Trie high_freq_kmer_trie;
	VectorBool ignore_kmer_bitmask;
	VectorBool ignore_window_bitmask;

	VectorBucketTables hash_tables;		// LSH min-hash index
} ref_t;

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
	//hash_t simhash;					// simhash / minhash (1bit)
	VectorMinHash minhashes;		// minhash vector
	char valid_minhash;

	// original mapping information
	int strand;
	uint32_t ref_pos_l;
	uint32_t ref_pos_r;
	
	// found ref match positions
	VectorU32 ref_bucket_id_matches_by_table;
	std::vector<VectorRefMatches> ref_matches;

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


void fasta2ref(const char *fastaFname, ref_t& ref);
void fastq2reads(const char *readsFname, reads_t& reads);
void print_read(read_t* read);
void parse_read_mapping(const char* read_name, unsigned int* ref_pos_l, unsigned int* ref_pos_r, int* strand);
void load_kmer_hist(const char* refFname, MapKmerCounts& hist, const uint32 max_count);
void store_kmer_hist(const char* refFname, const MapKmerCounts& hist);
void store_kmer_hist_stat(const char* refFname, const MapKmerCounts& hist);
void store_freq_kmers(const MapKmerCounts& hist);
void load_freq_kmers(const char* refFname, marisa::Trie& freq_trie, const uint32 max_count_threshold);

// index io
void store_ref_idx(const char* idxFname, const ref_t& ref);
void load_ref_idx(const char* idxFname, ref_t& ref);
void store_perm(const char* permFname, const VectorU32& perm);
void load_perm(const char* permFname, VectorU32& perm);
void store_hash_pads(const char* permFname, const VectorHash& perm);
void load_hash_pads(const char* permFname, VectorHash& perm);

// compression
#define CHARS_PER_SHORT 8   // number of chars in 16 bits
#define CHARS_PER_WORD 	16	// number of chars in 32 bits
#define BITS_PER_CHAR 	2
#define BITS_IN_SHORT 	16
#define BITS_IN_WORD 	32
#define BASE_IGNORE		4
#define KMER_HIST_SIZE16 (1ULL << 16) //65536
#define KMER_HIST_SIZE32 (1ULL << 32)

int pack_16(const char *seq, const int length, uint16_t *ret);
int pack_32(const char *seq, const int length, uint32_t *ret); 
void unpack_32(uint32 w, unsigned char *seq, const uint32 length);

#endif /*IO_H_*/
