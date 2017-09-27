#ifndef IO_H_
#define IO_H_

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
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

#define READ_BATCH_SIZE 100000
struct fastq_reader_t {
	seqan::SeqFileIn file_handle;
	int n_records;

	void open_file(const std::string& fname) {
		if (!seqan::open(file_handle, seqan::toCString(fname))) {
			std::cerr << "ERROR: Could not open FASTQ file: " << fname << "\n";
			exit(1);
		}
		n_records = 0;
	}
	
	// load FASTQ read records
	bool load_next_read(read_t& r) {
		if(!seqan::atEnd(file_handle)) {
			seqan::readRecord(r.name, r.seq, r.qual, file_handle);
			r.len = r.seq.size();
			for(uint32 i = 0; i < r.len; i++) {
				r.seq[i] = nt4_table[(int) r.seq[i]];
			}
			r.rc = r.seq;
			for(uint32 i = 0; i < r.len; i++) {
				r.rc[i] = nt4_complement[r.seq[r.len-i-1]];
			}
			r.rid = n_records;
			n_records++;
			return true;
		} else {
			return false;
		}
	}
	
	bool load_next_read_batch(reads_t& reads, const int read_batch_size) {
		int n_reads_loaded = 0;
		read_t r;
		while(n_reads_loaded < read_batch_size && this->load_next_read(r)) {
			reads.reads.push_back(r);
			n_reads_loaded++;
		}
		std::cout << "Loaded " << n_reads_loaded << " reads into current batch, total " << n_records << " reads loaded\n";
		return n_reads_loaded > 0;
	}
	
	void close_file() {
		seqan::close(file_handle);
	}
};

typedef enum {STORE, LOAD} precomp_contig_io_mode_t;
struct precomp_contig_io_t {
	precomp_contig_io_mode_t file_open_mode;
	std::ofstream file_handle_store;
	std::ifstream file_handle_load;
	
	void open_file(const char* fname, const precomp_contig_io_mode_t file_open_mode) {
		if(file_open_mode == STORE) {
			file_handle_store.open(fname, std::ios::out | std::ios::binary);
			if (!file_handle_store.is_open()) {
				printf("precomp_contig_io_t store: Cannot open file %s!\n", fname);
				exit(1);
			}
		} else { // LOAD
			file_handle_load.open(fname, std::ios::in | std::ios::binary);
			if (!file_handle_load.is_open()) {
				printf("precomp_contig_io_t load: Cannot open file %s!\n", fname);
				exit(1);
			}
		}
	}
	
	void store_precomp_contigs(reads_t& reads) {
		for(uint32 i = 0; i < reads.reads.size(); i++) {
			read_t* r = &reads.reads[i];
			if(!r->valid_minhash_f && !r->valid_minhash_rc) continue;
			uint32 ref_size = r->ref_matches.size();
			file_handle_store.write(reinterpret_cast<char*>(&ref_size), sizeof(r->ref_matches.size()));
			file_handle_store.write(reinterpret_cast<char*>(&(r->ref_matches[0])), r->ref_matches.size()*sizeof(ref_match_t));
		}
	}
	
	void load_precomp_contigs(reads_t& reads) {
		for(uint32 i = 0; i < reads.reads.size(); i++) {
			read_t* r = &reads.reads[i];
			if(!r->valid_minhash_f && !r->valid_minhash_rc) continue;
			uint32 ref_size;
			file_handle_load.read(reinterpret_cast<char*>(&ref_size), sizeof(r->ref_matches.size()));
			r->ref_matches.resize(ref_size);
			file_handle_load.read(reinterpret_cast<char*>(&(r->ref_matches[0])), r->ref_matches.size()*sizeof(ref_match_t));
			r->n_proc_contigs = ref_size;
			for(uint32 j = 0; j < r->ref_matches.size(); j++) {
				ref_match_t ref_contig = r->ref_matches[j];
				if(ref_contig.n_diff_bucket_hits > r->best_n_bucket_hits) {
					r->best_n_bucket_hits = ref_contig.n_diff_bucket_hits;
				}
			}
		}
	}
	
	void close_file() {
		if(file_handle_store.is_open()) file_handle_store.close();
		if(file_handle_load.is_open()) file_handle_load.close();
	}
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
void load_precomp_contigs(const char* fileName, reads_t& reads);
void store_precomp_contigs(const char* fileName, reads_t& reads);

// stats
void compute_and_store_kmer_hist32(const char* refFname, const char* seq, const seq_t seq_len, const index_params_t* params);
void compute_and_store_kmer_hist16(const char* refFname, const char* seq, const seq_t seq_len, const index_params_t* params);
void store_kmer_hist_stat(const char* refFname, const MapKmerCounts& hist);
bool load_freq_kmers(const char* refFname, std::set<uint64>& freq_kmers, const index_params_t* params);
bool load_freq_kmers(const char* refFname, VectorBool& freq_kmers_bitmap, MarisaTrie& freq_trie, const uint32 max_count_threshold);
void kmer_stats(const char* refFname);
void store_ref_index_stats(const char* refFname, const ref_t& ref, const index_params_t* params);
void ref_kmer_repeat_stats(const char* fastaFname, index_params_t* params, ref_t& ref);
void bin_repeat_stats(const char* fastaFname, index_params_t* params, ref_t& ref);

#endif /*IO_H_*/
