#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <iostream>
#include <fstream>

#include "io.h"
#include "hash.h"

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
				ref.seq.append(1, nt4_table[(int) c]);
			}
			c = (char) getc(fastaFile);
		}
	}
	printf("Done reading FASTA file. Total sequence length read = %zu\n", ref.seq.size());
	fclose(fastaFile);
}

// store the reference index
void store_ref_idx(const char* idxFname, ref_t& ref) {
//	std::ofstream file;
//	file.open(idxFname, std::ios::out | std::ios::app | std::ios::binary);
//	if (!file.is_open()) {
//		printf("store_ref_idx: Cannot open the IDX file %s!\n", idxFname);
//		exit(1);
//	}
//
//	// write the map
//	uint32 map_size = ref.windows_by_pos.size();
//	file.write(reinterpret_cast<char*>(&map_size), sizeof(map_size));
//	for (MapPos2Window::iterator it = ref.windows_by_pos.begin(); it != ref.windows_by_pos.end(); ++it) {
//		seq_t pos = it->first;
//		ref_win_t window = it->second;
//		uint32 n_minhashes = window.minhashes.size();
//
//	    file.write(reinterpret_cast<char*>(&pos), sizeof(pos));
//	    file.write(reinterpret_cast<char*>(&window.simhash), sizeof(window.simhash));
//	    file.write(reinterpret_cast<char*>(&n_minhashes), sizeof(n_minhashes));
//
//	    for(uint32 i = 0; i < n_minhashes; i++) {
//	    	hash_t minh = window.minhashes[i];
//	    	file.write(reinterpret_cast<char*>(&minh), sizeof(minh));
//	    }
//	}
//
//	// write the histogram
//	map_size = ref.kmer_hist.size();
//	file.write(reinterpret_cast<char*>(&map_size), sizeof(map_size));
//	for (MapKmerCounts::iterator it = ref.kmer_hist.begin(); it != ref.kmer_hist.end(); ++it) {
//		uint32 kmer = it->first;
//		seq_t count = it->second;
//		file.write(reinterpret_cast<char*>(&kmer), sizeof(kmer));
//		file.write(reinterpret_cast<char*>(&count), sizeof(count));
//	}
//	file.close();
}

// store the reference index
void load_ref_idx(const char* idxFname, ref_t& ref) {
//	std::ifstream file;
//	file.open(idxFname, std::ios::in | std::ios::binary);
//	if (!file.is_open()) {
//		printf("load_ref_idx: Cannot open the IDX file %s!\n", idxFname);
//		exit(1);
//	}
//
//	// read the windows map
//	uint32 map_size;
//	file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
//	for (uint32 i = 0; i < map_size; i++) {
//		seq_t pos;
//		ref_win_t window;
//		uint32 n_minhashes;
//		file.read(reinterpret_cast<char*>(&pos), sizeof(pos));
//		file.read(reinterpret_cast<char*>(&window.simhash), sizeof(window.simhash));
//		file.read(reinterpret_cast<char*>(&n_minhashes), sizeof(n_minhashes));
//
//		window.minhashes.resize(n_minhashes);
//		for(uint32 j = 0; j < n_minhashes; j++) {
//			hash_t minh;
//			file.read(reinterpret_cast<char*>(&minh), sizeof(minh));
//			window.minhashes[j] = minh;
//		}
//		ref.windows_by_pos.insert(std::pair<seq_t, ref_win_t>(pos, window));
//	}
//
//	// read the histogram
//	file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
//	for (uint32 i = 0; i < map_size; i++) {
//		uint32 kmer;
//		seq_t count;
//		file.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
//		file.read(reinterpret_cast<char*>(&count), sizeof(count));
//		ref.kmer_hist.insert(std::pair<uint32, seq_t> (kmer, count));
//	}
//	file.close();
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
