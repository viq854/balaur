#include <istream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <limits.h>
#include "io.h"
#include "types.h"
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

/* Reference I/O */

void fasta_error(const char* fastaFname) {
	printf("Error: File %s does not comply with the FASTA file format \n", fastaFname);
	exit(1);
}

/*void fasta2ref2(const char *fastaFname, ref_t& ref) {
	seqan::SeqFileIn file_handle;
	if (!seqan::open(file_handle, fastaFname)) {
		std::cerr << "ERROR: Could not open FASTA/FASTQ file: " << fastaFname << "\n";
		exit(1);
	}
	while(!seqan::atEnd(file_handle)) {
		std::string seq_name;
		std::string seq;
		ref.subsequence_offsets.push_back(ref.seq.size());
		seqan::readRecord(seq_name, seq, file_handle);
		ref.seq_names.push_back(seq_name);
		for(int i = 0; i < seq.size(); i++) {
			ref.seq.push_back(nt4_table[(int) seq[i]]);
		}
	}
	ref.len = ref.seq.size();
	seqan::close(file_handle);
	printf("Done reading FASTA file. Number of subsequences: %zu. Total sequence length read = %u\n", ref.subsequence_offsets.size(), ref.len);
}*/

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
		std::string seq_name;

		c = (char) getc(fastaFile);
		seq_name.push_back(c);		

		ref.subsequence_offsets.push_back(ref.seq.size());
		// sequence description line (> ...)
		while(c != '\n' && !feof(fastaFile)){
			c = (char) getc(fastaFile);
			seq_name.push_back(c);
		}
		seq_name.pop_back();
		if(feof(fastaFile)) fasta_error(fastaFname);
		ref.seq_names.push_back(seq_name);
		//std::cout << seq_name << "\n";

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

void compute_ref_repeat_mask(ref_t& ref) {
	const seq_t max_contig_len = params->max_matched_contig_len;///10;
	const seq_t n_contigs = ref.len - max_contig_len + 1;
	ref.contig_mask.resize(n_contigs);
	
	#pragma omp parallel for
	for(seq_t i = 0; i < n_contigs; i++) {
		for(int j = 0; j < max_contig_len; j++) {
			const uint16_t r = ref.precomputed_neighbor_repeats[i+j]; // distance to closest repeat
			const seq_t next_occ =  i + r;
			if(r == 0 || next_occ >= i + max_contig_len) continue; // unique kmer
			ref.contig_mask[i] = 1;
			break;
		}
	}
}

bool load_repeat_info(const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".rep.");
	fname += std::to_string(params->k2);
	fname += std::to_string(params->kmer_hashing_alg);
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("load_repeat_info: Could not open the file %s!\n", fname.c_str());
		return false;
	}
	ref.precomputed_neighbor_repeats.resize(ref.len - params->k2 + 1);
	file.read(reinterpret_cast<char*>(&ref.precomputed_neighbor_repeats[0]), (ref.len - params->k2 + 1)*sizeof(ref.precomputed_neighbor_repeats[0]));
	
	//const seq_t max_contig_len = params->max_matched_contig_len;
        //const seq_t n_contigs = ref.len - max_contig_len + 1;
	//ref.contig_mask.resize(n_contigs);
 	//file.read(reinterpret_cast<char*>(&ref.contig_mask[0]), ref.contig_mask.size()*sizeof(ref.contig_mask[0]));
	//file.close();
	return true;
}

/*bool load_repeat_local(const char* refFname, ref_t& ref, const index_params_t* params) {
        std::string fname(refFname);
        fname += std::string(".local_rep_map.");
        fname += std::to_string(params->k2);
        fname += std::to_string(params->kmer_hashing_alg);
        std::ifstream file;
        file.open(fname.c_str(), std::ios::in | std::ios::binary);
        if (!file.is_open()) {
                return false;
        }
	seq_t n_repeats;
	file.read(reinterpret_cast<char*>(&n_repeats), sizeof(n_repeats));
	
	std::cout << "N repeats: " << n_repeats << "\n";
	for(uint64_t i = 0; i < n_repeats; i++) {
		uint32 r;
        	file.read(reinterpret_cast<char*>(&r), sizeof(r));
		ref.repeats.insert((uint32) r);
		//ref.repeats_vec.push_back(r);
	}
	//std::sort(ref.repeats_vec.begin(), ref.repeats_vec.end());
	file.close();
        return true;
}*/

void compute_store_kmer2_hashes(const char* refFname, ref_t& ref, const index_params_t* params) {
	ref.precomputed_kmer2_hashes.resize(ref.len - params->k2 + 1);
	//#pragma omp parallel for
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
	compute_ref_repeat_mask(ref);

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
	file.write(reinterpret_cast<char*>(&ref.contig_mask[0]), ref.contig_mask.size()*sizeof(ref.contig_mask[0]));

	file.close();
}

/*void compute_store_repeat_local(const char* refFname, ref_t& ref, const index_params_t* params) {
        std::string fname(refFname);
        fname += std::string(".local_rep_map.");
        fname += std::to_string(params->k2);
        fname += std::to_string(params->kmer_hashing_alg);
        std::ofstream file;
        file.open(fname.c_str(), std::ios::out | std::ios::binary);
        if (!file.is_open()) {
                printf("compute_store_repeat_local: Cannot open the file %s!\n", fname.c_str());
                exit(1);
        }

	//#pragma omp parallel for
	seq_t n_repeats = 0;
        for (seq_t i = 0; i < ref.precomputed_kmer2_hashes.size(); i++) {
                uint32 k = (uint32) ref.precomputed_kmer2_hashes[i];
		if(ref.repeats.find(k) != ref.repeats.end()) continue;
                for(seq_t j = 1; j < params->ref_window_size; j++) {
                        if((i+j) >= ref.precomputed_kmer2_hashes.size()) break;
                        if(k == (uint32) ref.precomputed_kmer2_hashes[i + j]) {
                                ref.repeats.insert(k);
				n_repeats++;
				break;
                        }
                }
        }
	file.write(reinterpret_cast<char*>(&n_repeats), sizeof(n_repeats));
	for (auto& k: ref.repeats) {
		file.write(reinterpret_cast<const char*>(&k), sizeof(k));
	} 
        file.close();
}*/

bool load_kmer2_hashes(const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".hash.");
	fname += std::to_string(params->k2);
	fname += std::string(".alg.");
	fname += std::to_string(params->kmer_hashing_alg);
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
                printf("load_kmer2_hashes: Could not open the file %s!\n", fname.c_str());
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
		r.rid = reads.reads.size();
		reads.reads.push_back(r);
	}
	fclose(readsFile);
}

void store_precomp_contigs(const char* fileName, reads_t& reads) {
	std::string fname(fileName);
	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::app);
	if (!file.is_open()) {
		printf("store_or_load_contigs: Cannot open file %s!\n", fname.c_str());
		exit(1);
	}
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->valid_minhash_f && !r->valid_minhash_rc) continue;
		uint32 ref_size = r->ref_matches.size();
		file.write(reinterpret_cast<char*>(&ref_size), sizeof(r->ref_matches.size()));
		file.write(reinterpret_cast<char*>(&(r->ref_matches[0])), r->ref_matches.size()*sizeof(ref_match_t));
	}
}

struct ref_match_old_t {
	uint32_t pos;
	uint32 len;
	bool rc;
	int n_diff_bucket_hits;
};

void load_precomp_contigs(const char* fileName, reads_t& reads) {
	std::string fname(fileName);
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("store_or_load_contigs: Cannot open file %s!\n", fname.c_str());
		exit(1);
	}
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->valid_minhash_f && !r->valid_minhash_rc) continue;
		uint32 ref_size;
		file.read(reinterpret_cast<char*>(&ref_size), sizeof(r->ref_matches.size()));
		r->ref_matches.resize(ref_size);
		//std::vector<ref_match_old_t> matches(ref_size); 
		//file.read(reinterpret_cast<char*>(&(matches[0])), matches.size()*sizeof(ref_match_old_t));
		//for(int x = 0; x < ref_size; x++) {
		//	r->ref_matches[x] = ref_match_t(matches[x].pos, matches[x].len, matches[x].rc, matches[x].n_diff_bucket_hits);
		//}
		file.read(reinterpret_cast<char*>(&(r->ref_matches[0])), r->ref_matches.size()*sizeof(ref_match_t));
		r->n_proc_contigs = ref_size;
		for(uint32 j = 0; j < r->ref_matches.size(); j++) {
			ref_match_t ref_contig = r->ref_matches[j];
			if(ref_contig.n_diff_bucket_hits > r->best_n_bucket_hits) {
				r->best_n_bucket_hits = ref_contig.n_diff_bucket_hits;
			}
		}
	}
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
