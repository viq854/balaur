#ifndef IO_H_
#define IO_H_
#include "types.h"

#define INIT_REFLEN_ALLOC			262144
#define INIT_READLEN_ALLOC			100
#define MAX_SEQ_NAME_LEN 			256
#define MAX_NUM_READS 				1000000

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

typedef struct {
	seq_t pos;
	simhash_t simhash;
} ref_win_t;

typedef struct {
	seq_t len;
	char* seq;
	ref_win_t* windows;
	seq_t num_windows;
	int* hist;
} ref_t;

typedef struct {
	// read length
	int len;
	// read name
	char name[MAX_SEQ_NAME_LEN+1];
	// read sequence
	char* seq;
	// reverse complement sequence
	char* rc;
	// quality scores
	char* qual;
	
	// simhash fingerprint
	simhash_t simhash;
	int simhash_popc;

	// original mapping information
	int strand;
	uint64_t ref_pos_l;
	uint64_t ref_pos_r;
} read_t;

// collection of reads
typedef struct {
	// number of reads
	unsigned int count;
	read_t* reads;
	int* hist;
} reads_t;


ref_t* fasta2ref(char *fastaFname);
reads_t* fastq2reads(char *readsFname);

void print_read(read_t* read);
void free_reads(reads_t* reads);

void parse_read_mapping(char* read_name, unsigned int* ref_pos_l, unsigned int* ref_pos_r, int* strand); 

// index io
void store_ref_idx(ref_t* ref, const char* idxFname);
ref_t* load_ref_idx(const char* idxFname);

// compression
int pack_16(const char *seq, const int length, uint16_t* err);

#endif /*IO_H_*/
