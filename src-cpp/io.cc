#include <istream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <omp.h>
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

#define NUM_HIST_BUCKETS 1000000
void store_kmer_hist_stat(const char* refFname, const MapKmerCounts& hist) {
	std::string fname(refFname);
	fname += std::string(".kmer_hist_stats");
	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::app);
	if (!file.is_open()) {
		printf("store_kmer_hist: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}

	VectorU32 freq_buckets(NUM_HIST_BUCKETS);
	for (MapKmerCounts::const_iterator it = hist.begin(); it != hist.end(); ++it) {
		//uint32 kmer = it->first;
		seq_t count = it->second;
		if(count >= NUM_HIST_BUCKETS) {
			freq_buckets[NUM_HIST_BUCKETS-1]++;
		} else {
			freq_buckets[count]++;
		}
	}
	for(uint32 i = 0; i < NUM_HIST_BUCKETS; i++) {
		file << i << "\t" << freq_buckets[i] << "\n";
	}

	file.close();
}

void store_kmer_hist(const char* refFname, const MapKmerCounts& hist) {
	std::string fname(refFname);
	fname += std::string(".kmer_hist");

	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);

	if (!file.is_open()) {
		printf("store_kmer_hist: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}

	uint32 map_size = hist.size();
	file.write(reinterpret_cast<char*>(&map_size), sizeof(map_size));
	for (MapKmerCounts::const_iterator it = hist.begin(); it != hist.end(); ++it) {
		uint32 kmer = it->first;
		seq_t count = it->second;
		file.write(reinterpret_cast<char*>(&kmer), sizeof(kmer));
		file.write(reinterpret_cast<char*>(&count), sizeof(count));
	}
	file.close();
}

bool kmer_has_zero(uint32 kmer) {
	uint32 mask = (3 << (BITS_IN_WORD - BITS_PER_CHAR));
	for (int i = 0; i < 16; i++) {
		if (!(kmer & mask)) {
			return true;
		}
		mask >>= BITS_PER_CHAR;
	}
	return false;
}

void load_freq_kmers(const char* refFname, VectorBool& freq_kmers_bitmap, MarisaTrie& freq_trie, const uint32 max_count_threshold) {
	std::string fname(refFname);
	fname += std::string(".kmer_hist");

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);

	if (!file.is_open()) {
		printf("load_kmer_hist: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}

	freq_kmers_bitmap.resize(UINT_MAX + 1L);
	uint32 kmer, count;
	uint32 filtered = 0;
	uint32 tot_filtered = 0;
	int map_size;
	file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
#if(USE_MARISA)
	marisa::Keyset keys;
#endif
	while(map_size >= 0) {
		file.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
		file.read(reinterpret_cast<char*>(&count), sizeof(count));
		if(count >= max_count_threshold) {
			tot_filtered++;
			if (!kmer_has_zero(kmer)) {
				freq_kmers_bitmap[kmer] = true;
				filtered++;
			}
#if(USE_MARISA)
			unsigned char* seq = (unsigned char*) malloc(17*sizeof(char));
			unpack_32(kmer, seq, 16);
			seq[16] = '\0';
			keys.push_back((const char*) seq);
			free(seq);
#endif
		}
		map_size--;
	}
#if(USE_MARISA)
	freq_trie.build(keys, 0);
#endif
	file.close();
	printf("Filtered %u kmers, tot %u \n", filtered, tot_filtered);
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

void compute_store_kmer2_hashes(const char* refFname, ref_t& ref, const index_params_t* params) {
	ref.precomputed_kmer2_hashes.resize(ref.len - params->k2 + 1);
	#pragma omp parallel for
	for (seq_t pos = 0; pos < ref.len - params->k2 + 1; pos++) {
		ref.precomputed_kmer2_hashes[pos] = CityHash32(&ref.seq[pos], params->k2);
	}
	std::string fname(refFname);
	fname += std::string(".k2hash.");
	fname += std::to_string(params->k2);
	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);

	if (!file.is_open()) {
		printf("compute_store_kmer2_hashes: Cannot open the file %s!\n", fname.c_str());
		exit(1);
	}
	for (uint32 i = 0; i < ref.len - params->k2 + 1; i++) {
		file.write(reinterpret_cast<char*>(&ref.precomputed_kmer2_hashes[i]), sizeof(ref.precomputed_kmer2_hashes[i]));
	}
	file.close();
}

bool load_kmer2_hashes(const char* refFname, ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".k2hash.");
	fname += std::to_string(params->k2);
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		return false;
	}
	ref.precomputed_kmer2_hashes.resize(ref.len - params->k2 + 1);
	for(seq_t pos = 0; pos < ref.len - params->k2 + 1; pos++) {
		file.read(reinterpret_cast<char*>(&ref.precomputed_kmer2_hashes[pos]), sizeof(ref.precomputed_kmer2_hashes[pos]));
	}
	file.close();
	return true;
}

// store the reference index
void store_ref_idx(const char* refFname, const ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx.");
	fname += std::string("h");
	fname += std::to_string(params->h);
	fname += std::string("T");
	fname += std::to_string(params->n_tables);
	fname += std::string("b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("w");
	fname += std::to_string(params->ref_window_size);
	fname += std::string("__p");
        fname += std::to_string(params->n_buckets_pow2);

	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);
	if (!file.is_open()) {
		printf("store_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		exit(1);
	}

	uint32 size = ref.hash_tables.size();
	file.write(reinterpret_cast<char*>(&size), sizeof(size));
	for(uint32 i = 0; i < ref.hash_tables.size(); i++) {
		const buckets_t& buckets = ref.hash_tables[i];
		// indices
		file.write(reinterpret_cast<const char*>(&buckets.n_buckets), sizeof(buckets.n_buckets));
		for(uint32 j = 0; j < buckets.n_buckets; j++) {
			file.write(reinterpret_cast<const char*>(&buckets.bucket_indices[j]), sizeof(uint32));
		}
		// data
		size = buckets.next_free_bucket_index;
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		for(uint32 j = 0; j < buckets.n_buckets; j++) {
			if(buckets.bucket_indices[j] == buckets.n_buckets) {
				continue;
			}
			const VectorSeqPos& bucket = buckets.buckets_data_vectors[buckets.bucket_indices[j]];
			size = bucket.size();//buckets.bucket_sizes[buckets.bucket_indices[j]];
			file.write(reinterpret_cast<const char*>(&size), sizeof(size));
			for(uint32 k = 0; k < size; k++) {
				file.write(reinterpret_cast<const char*>(&bucket[k].pos), sizeof(seq_t));
				file.write(reinterpret_cast<const char*>(&bucket[k].len), sizeof(len_t));
			}
		}
	}
	file.close();
}

void store_ref_index_stats(const char* refFname, const ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx.");
	fname += std::string("h");
	fname += std::to_string(params->h);
	fname += std::string("T");
	fname += std::to_string(params->n_tables);
	fname += std::string("b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("w");
	fname += std::to_string(params->ref_window_size);
	fname += std::string("__p");
        fname += std::to_string(params->n_buckets_pow2);
	fname += std::string("__stats");

	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::app);
	if (!file.is_open()) {
		printf("store_ref_idx: Cannot open the STATS IDX file %s!\n", fname.c_str());
		exit(1);
	}

	for(uint32 i = 0; i < ref.hash_tables.size(); i++) {
		const buckets_t& buckets = ref.hash_tables[i];
		for(uint32 j = 0; j < buckets.n_buckets; j++) {
			if(buckets.bucket_indices[j] == buckets.n_buckets) {
				// table id, bucket id, size, average contig len
				file << i << "," << j << "," << "0" << "," << "0" << "\n";
				continue;
			}
			const VectorSeqPos& bucket = buckets.buckets_data_vectors[buckets.bucket_indices[j]];
			uint32 size = bucket.size();
			uint32 len_avg = 0;
			for(uint32 k = 0; k < size; k++) {
				len_avg += bucket[k].len;
			}
			if(bucket.size() > 0) {
				len_avg = len_avg/bucket.size();
			}
			// table id, bucket id, size, average contig len
			file << i << "," << j << "," << size << "," << len_avg << "\n";
		}
	}
	file.close();
}

// load the reference index buckets
void load_ref_idx(const char* refFname, ref_t& ref, index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx.");
	fname += std::string("h");
	fname += std::to_string(params->h);
	fname += std::string("T");
	fname += std::to_string(params->n_tables);
	fname += std::string("b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("w");
	fname += std::to_string(params->ref_window_size);
        fname += std::string("__p");
	fname += std::to_string(params->n_buckets_pow2);

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("load_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		std::cerr << "Error: " << strerror(errno);
		exit(1);
	}

	file.read(reinterpret_cast<char*>(&params->n_tables), sizeof(params->n_tables));
	ref.hash_tables.resize(params->n_tables);
	for(uint32 i = 0; i < params->n_tables; i++) {
		buckets_t* buckets = &ref.hash_tables[i];
		file.read(reinterpret_cast<char*>(&buckets->n_buckets), sizeof(buckets->n_buckets));
		buckets->bucket_indices.resize(buckets->n_buckets);
		for(uint32 j = 0; j < buckets->n_buckets; j++) {
			file.read(reinterpret_cast<char*>(&buckets->bucket_indices[j]), sizeof(buckets->bucket_indices[j]));
		}
		// note: no need to initialize the locks
		// data
		file.read(reinterpret_cast<char*>(&buckets->next_free_bucket_index), sizeof(buckets->next_free_bucket_index));
		buckets->buckets_data_vectors.resize(buckets->next_free_bucket_index);
		for(uint32 j = 0; j < buckets->n_buckets; j++) {
			if(buckets->bucket_indices[j] == buckets->n_buckets) {
				continue;
			}
			VectorSeqPos& bucket = buckets->buckets_data_vectors[buckets->bucket_indices[j]];
			uint32 size;
			file.read(reinterpret_cast<char*>(&size), sizeof(size));
			bucket.resize(size);
			// note: bucket size can now be the length of the vector
			for(uint32 k = 0; k < size; k++) {
				file.read(reinterpret_cast<char*>(&bucket[k].pos), sizeof(seq_t));
				file.read(reinterpret_cast<char*>(&bucket[k].len), sizeof(len_t));
			}
		}
	}
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

	for(uint32 i = 0; i < ref.hash_tables.size(); i++) {
		buckets_t& buckets = ref.hash_tables[i];
		for(uint32 j = 0; j < buckets.n_buckets; j++) {
			file.write(reinterpret_cast<const char*>(&buckets.per_thread_bucket_indices[tid][j]), sizeof(buckets.per_thread_bucket_indices[tid][j]));
		}
		for(uint32 j = 0; j < buckets.n_buckets; j++) {
			if(buckets.per_thread_bucket_indices[tid][j] == buckets.n_buckets) {
				continue;
			}
			VectorSeqPos& bucket = buckets.per_thread_buckets_data_vectors[tid][buckets.per_thread_bucket_indices[tid][j]];
			uint32 size = buckets.per_thread_bucket_sizes[tid][buckets.per_thread_bucket_indices[tid][j]];
			file.write(reinterpret_cast<const char*>(&size), sizeof(size));
			for(uint32 k = 0; k < size; k++) {
				file.write(reinterpret_cast<const char*>(&bucket[k].pos), sizeof(seq_t));
				file.write(reinterpret_cast<const char*>(&bucket[k].len), sizeof(len_t));
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
			buckets_t* buckets = &ref.hash_tables[i];
			for(uint32 j = 0; j < buckets->n_buckets; j++) {
				file.read(reinterpret_cast<char*>(&buckets->per_thread_bucket_indices[tid][j]), sizeof(buckets->per_thread_bucket_indices[tid][j]));
			}
			for(uint32 b = 0; b < buckets->n_buckets; b++) {
				if(buckets->per_thread_bucket_indices[tid][b] == buckets->n_buckets) {
					continue;
				}
				uint32 global_bucket_index = buckets->bucket_indices[b];
				if(global_bucket_index == buckets->n_buckets) {
					global_bucket_index = buckets->next_free_bucket_index;
					buckets->next_free_bucket_index++;
					buckets->bucket_indices[b] = global_bucket_index;
				}
				VectorSeqPos& global_bucket = buckets->buckets_data_vectors[global_bucket_index];
				uint32 size;
				file.read(reinterpret_cast<char*>(&size), sizeof(size));
				for(uint32 k = 0; k < size; k++) {
					loc_t w;
					file.read(reinterpret_cast<char*>(&w.pos), sizeof(seq_t));
					file.read(reinterpret_cast<char*>(&w.len), sizeof(len_t));
					global_bucket.push_back(w);
				}
			}
		}
	}
	file.close();
}

void store_perm( const char* permFname, const VectorU32& perm) {
	std::ofstream file;
	file.open(permFname, std::ios::out | std::ios::app | std::ios::binary);
	if (!file.is_open()) {
		printf("store_perm: Cannot open the perm file %s!\n", permFname);
		exit(1);
	}
	uint32 size = perm.size();
	file.write(reinterpret_cast<char*>(&size), sizeof(size));
	for (uint32 i = 0; i < size; i++) {
		uint32 p = perm[i];
		file.write(reinterpret_cast<char*>(&p), sizeof(p));
	}
	file.close();
}

void load_perm(const char* permFname, VectorU32& perm) {
	std::ifstream file;
	file.open(permFname, std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("store_perm: Cannot open the perm file %s!\n", permFname);
		exit(1);
	}
	uint32 size;
	file.read(reinterpret_cast<char*>(&size), sizeof(size));
	perm.resize(size);
	for (uint32 i = 0; i < size; i++) {
		uint32 p = perm[i];
		file.read(reinterpret_cast<char*>(&p), sizeof(p));
		perm[i] = p;
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
