#include <istream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <omp.h>
#include <limits.h>
#include "io.h"
#include "types.h"

/* Reference I/O */

void fasta_error(const char* fastaFname) {
	printf("Error: File %s does not comply with the FASTA file format \n", fastaFname);
	exit(1);
}

// reads the sequence data from the FASTA file
void fasta2ref(const char *fastaFname, ref_t& ref) {
	FILE* fastaFile = (FILE*) fopen(fastaFname, "r");
	if (fastaFile == NULL) {
		printf("fasta2ref: Cannot open FASTA file: %s!\n", fastaFname);
		exit(1);
	}
	char c = (char) getc(fastaFile);
	if(c != '>') fasta_error(fastaFname);
	while(!feof(fastaFile)) {
		c = (char) getc(fastaFile);

		ref.subsequence_offsets.push_back(ref.seq.size());
		// sequence description line (> ...)
		while(c != '\n' && !feof(fastaFile)){
			c = (char) getc(fastaFile);
		}
		if(feof(fastaFile)) fasta_error(fastaFname);

		// sequence data
		while(c != '>' && !feof(fastaFile)){
			if (c != '\n'){
				if (c >= 'a' && c <= 'z'){
					c += 'A'-'a';
				}
				ref.seq.append(1, nt4_table[(int) c]);
			}
			c = (char) getc(fastaFile);
		}
	}
	ref.len = ref.seq.size();
	printf("Done reading FASTA file. Number of subsequences: %zu. Total sequence length read = %u\n", ref.subsequence_offsets.size(), ref.len);
	fclose(fastaFile);
}

void store_valid_window_mask(const char* refFname, const ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".window_mask.");
	fname += std::to_string(params->ref_window_size);

	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);

	if (!file.is_open()) {
		printf("store_valid_window_mask: Cannot open the mask file %s!\n", fname.c_str());
		exit(1);
	}
	for (uint32 i = 0; i < ref.ignore_window_bitmask.size(); i++) {
		if(ref.ignore_window_bitmask[i]) {
			char b = '1';
			file.write(reinterpret_cast<char*>(&b), sizeof(char));
		} else {
			char b = '0';
			file.write(reinterpret_cast<char*>(&b), sizeof(char));
		}

	}
	file.close();
}

bool load_valid_window_mask(const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".window_mask.");
	fname += std::to_string(params->ref_window_size);

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);

	if (!file.is_open()) {
		printf("load_valid_window_mask: Could not open the mask file %s!\n", fname.c_str());
		return false;
	}
	char b;
	ref.ignore_window_bitmask.resize(ref.len - params->ref_window_size + 1);
	for(seq_t pos = 0; pos < ref.len - params->ref_window_size + 1; pos++) {
		file.read(reinterpret_cast<char*>(&b), sizeof(char));
		if(b == '1') {
			ref.ignore_window_bitmask[pos] = 1;
		}
	}
	file.close();
	return true;
}

bool load_repeat_info(const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".rep.");
	fname += std::to_string(params->k2);
	fname += std::to_string(params->kmer_hashing_alg);
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		return false;
	}
	ref.precomputed_neighbor_repeats.resize(ref.len - params->k2 + 1);
	file.read(reinterpret_cast<char*>(&ref.precomputed_neighbor_repeats[0]), (ref.len - params->k2 + 1)*sizeof(ref.precomputed_neighbor_repeats[0]));
	file.close();
	return true;
}

void compute_store_kmer2_hashes(const char* refFname, ref_t& ref, const index_params_t* params) {
	ref.precomputed_kmer2_hashes.resize(ref.len - params->k2 + 1);
	#pragma omp parallel for
	for (seq_t pos = 0; pos < ref.len - params->k2 + 1; pos++) {
		switch(params->kmer_hashing_alg) {
			case SHA1_E:
				uint32_t hash[5];
				sha1_hash(reinterpret_cast<const uint8_t*>(&ref.seq[pos]), params->k2, hash);
				ref.precomputed_kmer2_hashes[pos] = ((uint64) hash[0] << 32 | hash[1]);
				break;
			case CITY_HASH64:
				ref.precomputed_kmer2_hashes[pos] = CityHash64(&ref.seq[pos], params->k2);
				break;
			case PACK64:
				pack_64(&ref.seq[pos], params->k2, &ref.precomputed_kmer2_hashes[pos]);
				break;
		}
	}
	std::string fname(refFname);
	fname += std::string(".hash.");
	fname += std::to_string(params->k2);
	fname += std::string(".alg.");
	fname += std::to_string(params->kmer_hashing_alg);
	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);
	if (!file.is_open()) {
		printf("compute_store_k2_hashes: Cannot open the file %s!\n", fname.c_str());
		exit(1);
	}
	for (uint32 i = 0; i < ref.len - params->k2 + 1; i++) {
		file.write(reinterpret_cast<char*>(&ref.precomputed_kmer2_hashes[i]), sizeof(ref.precomputed_kmer2_hashes[i]));
	}
	file.close();
}

void compute_store_repeat_info(const char* refFname, ref_t& ref, const index_params_t* params) {
	ref.precomputed_neighbor_repeats.resize(ref.len - params->k2 + 1);
	#pragma omp parallel for
	for (seq_t i = 0; i < ref.len - params->k2 + 1; i++) {
		uint64 k = ref.precomputed_kmer2_hashes[i];
		for(seq_t j = 1; j < MAX_LOC_LEN; j++) {
			if((i+j) == ref.precomputed_kmer2_hashes.size()) break;
			if(k == ref.precomputed_kmer2_hashes[i + j]) {
				ref.precomputed_neighbor_repeats[i] = j;
				break;
			}
		}
	}
	std::string fname(refFname);
	fname += std::string(".rep.");
	fname += std::to_string(params->k2);
	fname += std::to_string(params->kmer_hashing_alg);
	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);
	if (!file.is_open()) {
		printf("compute_store_repeat_info: Cannot open the file %s!\n", fname.c_str());
		exit(1);
	}
	for (uint32 i = 0; i < ref.len - params->k2 + 1; i++) {
		file.write(reinterpret_cast<char*>(&ref.precomputed_neighbor_repeats[i]), sizeof(ref.precomputed_neighbor_repeats[i]));
	}
	file.close();
}

bool load_kmer2_hashes(const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".hash.");
	fname += std::to_string(params->k2);
	fname += std::string(".alg.");
	fname += std::to_string(params->kmer_hashing_alg);
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		return false;
	}
	ref.precomputed_kmer2_hashes.resize(ref.len - params->k2 + 1);
	file.read(reinterpret_cast<char*>(&ref.precomputed_kmer2_hashes[0]), (ref.len - params->k2 + 1)*sizeof(ref.precomputed_kmer2_hashes[0]));
	file.close();
	return true;
}

void store_ref_idx_flat(const char* refFname, const ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx_flat.");
	fname += std::string("h");
	fname += std::to_string(params->h);
	fname += std::string("_T");
	fname += std::to_string(params->n_tables);
	fname += std::string("_b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("_w");
	fname += std::to_string(params->ref_window_size);
	fname += std::string("_p");
    fname += std::to_string(params->n_buckets_pow2);
    fname += std::string("_k");
    fname += std::to_string(params->k);
    fname += std::string("_H");
    fname += std::to_string(params->max_count);

	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);
	if (!file.is_open()) {
		printf("store_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		exit(1);
	}

	uint64 total_num_entries = ref.index.buckets_data.size();
	file.write(reinterpret_cast<char*>(&total_num_entries), sizeof(total_num_entries));
	file.write(reinterpret_cast<const char*>(&ref.index.buckets_data[0]), total_num_entries*sizeof(loc_t));
	file.write(reinterpret_cast<const char*>(&ref.index.bucket_offsets[0]), (params->n_tables*params->n_buckets+1)*sizeof(uint64));
	file.close();
}

void load_ref_idx_flat(const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx_flat.");
	fname += std::string("h");
	fname += std::to_string(params->h);
	fname += std::string("_T");
	fname += std::to_string(params->n_tables);
	fname += std::string("_b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("_w");
	fname += std::to_string(params->ref_window_size);
	fname += std::string("_p");
	fname += std::to_string(params->n_buckets_pow2);
	fname += std::string("_k");
	fname += std::to_string(params->k);
	fname += std::string("_H");
	fname += std::to_string(params->max_count);

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("load_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		std::cerr << "Error: " << strerror(errno);
		exit(1);
	}

	uint64 total_num_bucket_entries;
	file.read(reinterpret_cast<char*>(&total_num_bucket_entries), sizeof(total_num_bucket_entries));
	ref.index.bucket_offsets.resize(params->n_tables*params->n_buckets+1);
	ref.index.buckets_data.resize(total_num_bucket_entries);
	file.read(reinterpret_cast<char*>(&ref.index.buckets_data[0]), total_num_bucket_entries*sizeof(loc_t));
	file.read(reinterpret_cast<char*>(&ref.index.bucket_offsets[0]), (params->n_tables*params->n_buckets+1)*sizeof(uint64));
	file.close();
}

// store the reference index
void store_ref_idx(const char* refFname, const ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx.");
	fname += std::string("h");
	fname += std::to_string(params->h);
	fname += std::string("_T");
	fname += std::to_string(params->n_tables);
	fname += std::string("_b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("_w");
	fname += std::to_string(params->ref_window_size);
	fname += std::string("_p");
    fname += std::to_string(params->n_buckets_pow2);
    fname += std::string("_k");
    fname += std::to_string(params->k);
    fname += std::string("_H");
    fname += std::to_string(params->max_count);

	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);
	if (!file.is_open()) {
		printf("store_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		exit(1);
	}

	uint64 total_num_entries = 0;
	for(uint32 i = 0; i < params->n_tables; i++) {
		total_num_entries += ref.mutable_index.per_table_buckets[i].n_entries;
	}
	file.write(reinterpret_cast<char*>(&total_num_entries), sizeof(total_num_entries));
	for(uint32 i = 0; i < params->n_tables; i++) {
		const buckets_t& buckets = ref.mutable_index.per_table_buckets[i];
		for(uint32 j = 0; j < params->n_buckets; j++) {
			const VectorSeqPos& bucket = buckets.buckets_data_vectors[j];
			uint32 size = bucket.size();
			file.write(reinterpret_cast<const char*>(&size), sizeof(size));
			file.write(reinterpret_cast<const char*>(&bucket[0]), size*sizeof(loc_t));
		}
	}
	file.close();
}

// load the reference index buckets
void load_ref_idx(const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx.");
	fname += std::string("h");
	fname += std::to_string(params->h);
	fname += std::string("_T");
	fname += std::to_string(params->n_tables);
	fname += std::string("_b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("_w");
	fname += std::to_string(params->ref_window_size);
	fname += std::string("_p");
	fname += std::to_string(params->n_buckets_pow2);
	fname += std::string("_k");
	fname += std::to_string(params->k);
	fname += std::string("_H");
	fname += std::to_string(params->max_count);

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("load_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		std::cerr << "Error: " << strerror(errno);
		exit(1);
	}

	uint64 total_num_bucket_entries;
	file.read(reinterpret_cast<char*>(&total_num_bucket_entries), sizeof(total_num_bucket_entries));
	ref.index.bucket_offsets.resize(params->n_tables*params->n_buckets+1);
	ref.index.buckets_data.resize(total_num_bucket_entries);
	std::cout << "Total number of contig entries in the index: " << total_num_bucket_entries << "\n";

	uint64 bucket_idx = 0;
	for(uint32 i = 0; i < params->n_tables; i++) {
		for(uint32 j = 0; j < params->n_buckets; j++) {
			ref.index.bucket_offsets[i*params->n_buckets + j] = bucket_idx;
			uint32 size;
			file.read(reinterpret_cast<char*>(&size), sizeof(size));
			file.read(reinterpret_cast<char*>(&ref.index.buckets_data[bucket_idx]), size*sizeof(loc_t));
			std::sort(ref.index.buckets_data.begin() + bucket_idx, ref.index.buckets_data.begin() + bucket_idx + size, comp_loc());
			bucket_idx += size;
		}
	}
	ref.index.bucket_offsets[ref.index.bucket_offsets.size()-1] = bucket_idx;
	file.close();
}

void store_ref_idx_per_thread(const int tid, const bool first_entry, const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx_tid");
	fname += std::to_string(tid);

	std::ofstream file;
	if(first_entry) {
		file.open(fname.c_str(), std::ios::out | std::ios::binary);
	} else {
		file.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::app);
	}
	if (!file.is_open()) {
		printf("store_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		exit(1);
	}

	for(uint32 i = 0; i < params->n_tables; i++) {
		buckets_t& buckets = ref.mutable_index.per_table_buckets[i];
		for(uint32 j = 0; j < params->n_buckets; j++) {
			VectorSeqPos& bucket = buckets.per_thread_buckets_data_vectors[tid][j];
			uint32 size = buckets.per_thread_bucket_sizes[tid][j];
			file.write(reinterpret_cast<const char*>(&size), sizeof(size));
			for(uint32 k = 0; k < size; k++) {
				file.write(reinterpret_cast<const char*>(&bucket[k].pos), sizeof(seq_t));
				file.write(reinterpret_cast<const char*>(&bucket[k].len), sizeof(len_t));
				file.write(reinterpret_cast<const char*>(&bucket[k].hash), sizeof(minhash_t));
			}
			bucket.resize(0);
			bucket.shrink_to_fit();
		}
		std::fill(buckets.per_thread_bucket_sizes[tid].begin(), buckets.per_thread_bucket_sizes[tid].end(), 0);
	}
	file.close();
}

void load_ref_idx_per_thread(const int tid, const int nloads, const char* refFname, ref_t& ref, index_params_t* params) {
	if(nloads == 0) return;

	std::string fname(refFname);
	fname += std::string(".idx_tid");
	fname += std::to_string(tid);

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("load_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		std::cerr << "Error: " << strerror(errno) << "\n";
		exit(1);
	}

	for(int l = 0; l < nloads; l++) {
		for(uint32 i = 0; i < params->n_tables; i++) {
			buckets_t* buckets = &ref.mutable_index.per_table_buckets[i];
			for(uint32 b = 0; b < params->n_buckets; b++) {
				VectorSeqPos& global_bucket = buckets->buckets_data_vectors[b];
				uint32 size;
				file.read(reinterpret_cast<char*>(&size), sizeof(size));
				for(uint32 k = 0; k < size; k++) {
					loc_t w;
					file.read(reinterpret_cast<char*>(&w.pos), sizeof(seq_t));
					file.read(reinterpret_cast<char*>(&w.len), sizeof(len_t));
					file.read(reinterpret_cast<char*>(&w.hash), sizeof(minhash_t));
					global_bucket.push_back(w);
				}
			}
		}
	}
	file.close();
}

/* Reads I/O */

void fastq_error(const char* fastqFname) {
	printf("Error: File %s does not comply with the FASTQ file format \n", fastqFname);
	exit(1);
}

// loads the read sequences from the FASTQ file
void fastq2reads(const char *readsFname, reads_t& reads) {
	FILE *readsFile = (FILE*) fopen(readsFname, "r");
	if (readsFile == NULL) {
		printf("load_reads_fastq: Cannot open reads file: %s !\n", readsFname);
		exit(1);
	}

	reads.fname = readsFname;
	char c;
	while(!feof(readsFile)) {
		read_t r;
		c = (char) getc(readsFile);
		while(c != '@' && !feof(readsFile)) {
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) break;

		// line 1 (@ ...)
		c = (char) getc(readsFile);
		while(c != '\n' && !feof(readsFile)){
			r.name.append(1, c);
			c = (char) getc(readsFile);
		}
		r.name.append(1, '\0');
		if(feof(readsFile)) fastq_error(readsFname);

		while (c != '\n' && !feof(readsFile)) {
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) fastq_error(readsFname);

		// line 2 (sequence letters)
		c = (char) getc(readsFile);
		while (c != '\n' && !feof(readsFile)) {
			r.seq.append(1, nt4_table[(int) c]);
			c = (char) getc(readsFile);
		}
		r.len = r.seq.size();
		if(feof(readsFile)) fastq_error(readsFname);

		while (c != '+' && !feof(readsFile)) {
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) fastq_error(readsFname);

		// line 3 (+ ...)
		while(c != '\n' && !feof(readsFile)){
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) fastq_error(readsFname);

		// line 4 (quality values)
		uint32 qualLen = 0;
		c = (char) getc(readsFile);
		while(c != '\n' && !feof(readsFile)) {
			qualLen++;
			c = (char) getc(readsFile);
		}
		if(qualLen != r.len) {
			printf("Error: The number of quality score symbols does not match the length of the read sequence.\n");
			exit(1);
		}

		// compute the reverse complement
		for(uint32 i = 0; i < r.len; i++) {
			r.rc.append(1, nt4_complement[(int)r.seq.at(r.len-i-1)]);
		}
		reads.reads.push_back(r);
	}
	fclose(readsFile);
}

// assumes that reads were generated with wgsim
void parse_read_mapping(const char* read_name, unsigned int* seq_id, unsigned int* ref_pos_l, unsigned int* ref_pos_r, int* strand) {
    std::istringstream is((std::string(read_name)));
    std::string seqid;
    std::string refl;    
    std::string refr;
    std::getline(is, seqid, '_');
    std::getline(is, refl, '_');
    std::getline(is, refr, '_');
    if(seqid.compare("X") == 0) {
    	*seq_id = 23;
    } else {
    	*seq_id = atoi(seqid.c_str());
    }
    *ref_pos_l = atoi(refl.c_str());
    *ref_pos_r = atoi(refr.c_str());
}

void get_sim_read_info(const ref_t& ref, reads_t& reads) {
#if(SIM_EVAL)
	#pragma omp parallel for
	for (uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		parse_read_mapping(r->name.c_str(), &r->seq_id, &r->ref_pos_l, &r->ref_pos_r, &r->strand);
		r->seq_id = r->seq_id - 1;
		if(ref.subsequence_offsets.size() > 1) {
			r->ref_pos_l += ref.subsequence_offsets[r->seq_id]; // convert to global id
			r->ref_pos_r += ref.subsequence_offsets[r->seq_id];
		}
	}
#endif
}

void print_read(read_t* read) {
	printf("%s \n", read->name.c_str());
	for(uint32 i = 0; i < read->len; i++) {
		printf("%c", iupacChar[(int) read->seq[i]]);
	}
	printf("\n");
	for(uint32 i = 0; i < read->len; i++) {
		printf("%c", iupacChar[(int) read->rc[i]]);
	}
	printf("\n");
}

// --- Compression ---

// compress the seq of given length into 16 bits (using 2 bits per char)
// returns -1 if the seq contains bases that should be ignored (e.g. N or $)
int pack_16(const char *seq, const int length, uint16_t* ret) {
	uint16_t c = 0;
	for (int k = 0; k < length; k++) {
		if(seq[k] == BASE_IGNORE) {
			return -1;
		}
		c = c | (seq[k] << (BITS_IN_SHORT - (k+1) * BITS_PER_CHAR));
	}
	
	*ret = c;
	return 0;
}

void unpack_32(uint32 w, unsigned char *seq, const uint32 length) {
	for (uint32 k = 0; k < length; k++) {
		seq[k] = w >> (BITS_IN_WORD - BITS_PER_CHAR);
		w <<= BITS_PER_CHAR;
	}
}

int pack_32(const char *seq, const int length, uint32_t *ret) {
	uint32_t c = 0;
	for (int k = 0; k < length; k++) {
		if(seq[k] == BASE_IGNORE) {
			return -1;
		}
		c = c | (seq[k] << (BITS_IN_WORD - (k+1) * BITS_PER_CHAR));
	}
	*ret = c;
	return 0;
}


int pack_64(const char *seq, const int length, uint64 *ret) {
	uint64_t c = 0;
	for (int k = 0; k < length; k++) {
		if(seq[k] == BASE_IGNORE) {
			return -1;
		}
		c = c | ((seq[k] & 3ULL) << (BITS_IN_LWORD - (k+1) * BITS_PER_CHAR));
	}
	*ret = c;
	return 0;
}
