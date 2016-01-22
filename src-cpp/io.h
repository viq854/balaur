#ifndef IO_H_
#define IO_H_

#include "types.h"
#include "index.h"

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


void fasta2ref(const char *fastaFname, ref_t& ref);
void fastq2reads(const char *readsFname, reads_t& reads);
void print_read(read_t* read);
void parse_read_mapping(const char* read_name, unsigned int* seq_id, unsigned int* ref_pos_l, unsigned int* ref_pos_r, int* strand);
void get_sim_read_info(const ref_t& ref, reads_t& reads);
void store_valid_window_mask(const char* refFname, const ref_t& ref, const index_params_t* params);
bool load_valid_window_mask(const char* refFname, ref_t& ref, const index_params_t* params);

// index io
void store_ref_idx_flat(const char* refFname, const ref_t& ref, const index_params_t* params);
void load_ref_idx_flat(const char* refFname, ref_t& ref, const index_params_t* params);  
void store_ref_idx(const char* idxFname, const ref_t& ref, const index_params_t* params);
void load_ref_idx(const char* idxFname, ref_t& ref, const index_params_t* params);
void store_ref_idx_per_thread(const int tid, const bool first_entry, const char* refFname, ref_t& ref, const index_params_t* params);
void load_ref_idx_per_thread(const int tid, const int nloads, const char* refFname, ref_t& ref, index_params_t* params);
void compute_store_kmer2_hashes(const char* refFname, ref_t& ref, const index_params_t* params);
bool load_repeat_info(const char* refFname, ref_t& ref, const index_params_t* params);
void compute_store_repeat_info(const char* refFname, ref_t& ref, const index_params_t* params);
bool load_kmer2_hashes(const char* refFname, ref_t& ref, const index_params_t* params);
void mark_windows_to_discard(ref_t& ref, const index_params_t* params);
void mark_freq_kmers(ref_t& ref, const index_params_t* params);
void compute_store_repeat_local(const char* refFname, ref_t& ref, const index_params_t* params);
bool load_repeat_local(const char* refFname, ref_t& ref, const index_params_t* params);


// stats
void compute_and_store_kmer_hist32(const char* refFname, const char* seq, const seq_t seq_len, const index_params_t* params);
void compute_and_store_kmer_hist16(const char* refFname, const char* seq, const seq_t seq_len, const index_params_t* params);
void store_kmer_hist_stat(const char* refFname, const MapKmerCounts& hist);
void load_freq_kmers(const char* refFname, std::set<uint64>& freq_kmers, const index_params_t* params);
void load_freq_kmers(const char* refFname, VectorBool& freq_kmers_bitmap, MarisaTrie& freq_trie, const uint32 max_count_threshold);
void kmer_stats(const char* refFname);
void store_ref_index_stats(const char* refFname, const ref_t& ref, const index_params_t* params);
void ref_kmer_repeat_stats(const char* fastaFname, index_params_t* params, ref_t& ref);

// compression
#define CHARS_PER_SHORT 8   // number of chars in 16 bits
#define CHARS_PER_WORD 	16	// number of chars in 32 bits
#define BITS_PER_CHAR 	2
#define BITS_IN_SHORT 	16
#define BITS_IN_WORD 	32
#define BITS_IN_LWORD 	64
#define BASE_IGNORE		4
#define KMER_HIST_SIZE16 (1ULL << 16) //65536
#define KMER_HIST_SIZE32 (1ULL << 32)

int pack_16(const char *seq, const int length, uint16_t *ret);
int pack_32(const char *seq, const int length, uint32_t *ret); 
int pack_64(const char *seq, const int length, uint64 *ret);
void unpack_32(uint32 w, unsigned char *seq, const uint32 length);

#endif /*IO_H_*/
