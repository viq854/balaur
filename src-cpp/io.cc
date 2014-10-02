
#include <stdio.h>
#include <string.h>
#include <fstream>

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
				//if(ref.seq.size() < 100000000) {
					ref.seq.append(1, nt4_table[(int) c]);
				//}
			}
			c = (char) getc(fastaFile);
		}
	}
	ref.len = ref.seq.size();
	printf("Done reading FASTA file. Total sequence length read = %u\n", ref.len);
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
	file.open(fname.c_str(), std::ios::out | std::ios::app | std::ios::binary);

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

void load_freq_kmers(const char* refFname, marisa::Trie& freq_trie, const uint32 max_count_threshold) {
	std::string fname(refFname);
	fname += std::string(".kmer_hist");

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);

	if (!file.is_open()) {
		printf("load_kmer_hist: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}

	marisa::Keyset keys;
	uint32 kmer, count;
	uint32 filtered = 0;
	int map_size;
	file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
	while(map_size >= 0) {
		file.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
		file.read(reinterpret_cast<char*>(&count), sizeof(count));
		if(count >= max_count_threshold) {
			unsigned char* seq = (unsigned char*) malloc(17*sizeof(char));
			unpack_32(kmer, seq, 16);
			seq[16] = '\0';
			keys.push_back((const char*) seq);
			filtered++;
		}
		map_size--;
	}
	freq_trie.build(keys, 0);
	file.close();
	printf("Filtered %u kmers \n", filtered);
}


// store the reference index
void store_ref_idx(const char* refFname, const ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx");

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
			size = buckets.bucket_sizes[buckets.bucket_indices[j]];
			file.write(reinterpret_cast<const char*>(&size), sizeof(size));
			for(uint32 k = 0; k < size; k++) {
				file.write(reinterpret_cast<const char*>(&bucket[k]), sizeof(seq_t));
			}
		}
	}
	file.close();
}

// load the reference index buckets
void load_ref_idx(const char* refFname, ref_t& ref, index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx");

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("load_ref_idx: Cannot open the IDX file %s!\n", fname.c_str());
		cerr << "Error: " << strerror(errno);
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
				file.read(reinterpret_cast<char*>(&bucket[k]), sizeof(bucket[k]));
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

		// TODO: compute the reverse complement

		r.acc = 0;
		reads.reads.push_back(r);

		//if(reads.reads.size() > 1000) break;
	}
	fclose(readsFile);
}

// assumes that reads were generated with wgsim
void parse_read_mapping(const char* read_name, unsigned int* ref_pos_l, unsigned int* ref_pos_r, int* strand) {
    int token_index = 0;
    const char delimiters[] = "_";
    const char* read_name_dup = std::string(read_name).c_str();
    char* ptr;
    char* token = strtok_r((char *) read_name_dup, delimiters, &ptr);
    while (token != NULL) {
        if(token_index == 1) {
                sscanf(token, "%u", ref_pos_l);
        } else if(token_index == 2) {
                sscanf(token, "%u", ref_pos_r);
        } else if(token_index == 3) {
                *strand = (strcmp(token, "nm") == 0) ? 1 : 0;
        } else if(token_index > 3){
        	break;
        }
        token = strtok_r(NULL, delimiters, &ptr);
        token_index++;
    }
}

void print_read(read_t* read) {
	printf("%s \n", read->name.c_str());
	for(uint32 i = 0; i < read->len; i++) {
		printf("%c", iupacChar[(int) read->seq[i]]);
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
