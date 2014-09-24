#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include "types.h"

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
	VectorU32 bucket_sizes;
	VectorBuckets buckets_data_vectors;
};
typedef std::vector<buckets_t> VectorBucketTables;


// reference genome
typedef struct {
	std::string seq; 	// reference sequence data
	seq_t len;			// reference sequence length
	MapKmerCounts kmer_hist;				// kmer occurrence histogram
	MapKmerCounts high_freq_kmer_hist;

	// LSH index
	VectorBucketTables hash_tables;
} ref_t;

typedef struct {
	uint32_t len; 					// read length
	std::string name; 				// read name
	std::string seq;				// read sequence
	std::string rc;					// reverse complement sequence
	std::string qual;				// quality scores
	
	// LSH fingerprints
	hash_t simhash;					// simhash / minhash (1bit)
	VectorMinHash minhashes;		// minhash vector

	// original mapping information
	int strand;
	uint64_t ref_pos_l;
	uint64_t ref_pos_r;
	
	// found ref match positions
	VectorSeqPos ref_matches;
	MapPos2MinCount matched_window_counts;

	char acc; // DEBUG: whether read matched accurately
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

// index io
void store_ref_idx(const char* idxFname, ref_t& ref);
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

#endif /*IO_H_*/
