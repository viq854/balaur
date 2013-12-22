#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "io.h"
#include "hash.h"

/* Reference I/O */

void fasta_error(char* fastaFname) {
	printf("Error: File %s does not comply with the FASTA file format \n", fastaFname);
	exit(1);
}

// reads the sequence data from the FASTA file
ref_t* fasta2ref(char *fastaFname) {
	FILE * fastaFile = (FILE*) fopen(fastaFname, "r");
	if (fastaFile == NULL) {
		printf("fasta2ref: Cannot open FASTA file: %s!\n", fastaFname);
		exit(1);
	}

	ref_t* ref = (ref_t*) calloc(1, sizeof(ref_t));
	seq_t allocatedSeqLen = INIT_REFLEN_ALLOC;
	ref->len = 0;
	ref->seq = (char*) malloc(allocatedSeqLen * sizeof(char));

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
				// reallocate twice as much memory
				if (ref->len >= allocatedSeqLen) {
					allocatedSeqLen <<= 1;
					ref->seq = (char*)realloc(ref->seq, sizeof(char)*allocatedSeqLen);
				}
				*(ref->seq + ref->len) = nt4_table[(int) c];
				ref->len++;
			}
			c = (char) getc(fastaFile);
		}
		// add $ as a separator between chromosomes
		// to disallow spanning the boundary
		*(ref->seq + ref->len) = nt4_table[(int) '$'];
		ref->len++;
	}
	ref->len--; // to ignore the last $
	printf("Done reading FASTA file. Total sequence length read = %llu\n", ref->len);
	fclose(fastaFile);

	return ref;
}

void store_ref_idx(ref_t* ref, const char* idxFname) {
	FILE* idxFile = (FILE*) fopen(idxFname, "wb");
	if (idxFile == NULL) {
		printf("store_ref_idx: Cannot open the IDX file %s!\n", idxFname);
		exit(1);
	}
	fwrite(&ref->len, sizeof(seq_t), 1, idxFile);
	fwrite(&ref->num_windows, sizeof(seq_t), 1, idxFile);
	fwrite(&ref->hist_size, sizeof(seq_t), 1, idxFile);
	fwrite(ref->windows, sizeof(ref_win_t), ref->num_windows, idxFile);
	fwrite(ref->hist, sizeof(int), ref->hist_size, idxFile);
	fclose(idxFile);
}

ref_t* load_ref_idx(const char* idxFname) {
	FILE* idxFile = (FILE*) fopen(idxFname, "rb");
	if (idxFile == NULL) {
		printf("load_ref_idx: Cannot open the IDX file %s!\n", idxFname);
		exit(1);
	}
	ref_t* ref = (ref_t*) calloc(1, sizeof(ref_t));
	fread(&ref->len, sizeof(seq_t), 1, idxFile);
	fread(&ref->num_windows, sizeof(seq_t), 1, idxFile);
	fread(&ref->hist_size, sizeof(seq_t), 1, idxFile);
	ref->windows = (ref_win_t*) calloc(ref->num_windows, sizeof(ref_win_t));
	ref->hist = (int*) calloc(ref->hist_size, sizeof(int));
	if((ref->windows == 0) || (ref->hist == 0)) {
		printf("Could not allocate memory for the ref index. \n");
		exit(1);
	}
	fread(ref->windows, sizeof(ref_win_t), ref->num_windows, idxFile);
	fread(ref->hist, sizeof(int), ref->hist_size, idxFile);
	fclose(idxFile);
	
	return ref;
}

/* Reads I/O */

void fastq_error(char* fastqFname) {
	printf("Error: File %s does not comply with the FASTQ file format \n", fastqFname);
	exit(1);
}

// loads the read sequences from the FASTQ file
reads_t* fastq2reads(char *readsFname) {
	FILE *readsFile = (FILE*) fopen(readsFname, "r");
	if (readsFile == NULL) {
		printf("load_reads_fastq: Cannot open reads file: %s !\n", readsFname);
		exit(1);
	}
	reads_t *reads = (reads_t*) calloc(1, sizeof(reads_t));
	reads->reads = (read_t*) malloc(MAX_NUM_READS*sizeof(read_t));
	reads->count = 0;

	char c;
	while(!feof(readsFile)) {
		read_t* read = &(reads->reads[reads->count]);

		c = (char) getc(readsFile);
		while(c != '@' && !feof(readsFile)) {
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) break;

		// line 1 (@ ...)
		int seqNameLen = 0;
		c = (char) getc(readsFile);
		while(c != '\n' && seqNameLen < MAX_SEQ_NAME_LEN && !feof(readsFile)){
			read->name[seqNameLen] = c;
			seqNameLen++;
			c = (char) getc(readsFile);
		}
		read->name[seqNameLen] = '\0';
		if(feof(readsFile)) fastq_error(readsFname);

		while (c != '\n' && !feof(readsFile)) {
			c = (char) getc(readsFile);
		}
		if(feof(readsFile)) fastq_error(readsFname);

		// line 2 (sequence letters)
		int readLen = 0;
		int allocatedReadLen = INIT_READLEN_ALLOC;
		read->seq = (char*) malloc(allocatedReadLen*sizeof(char));
		read->rc = (char*) malloc(allocatedReadLen*sizeof(char));
		//read->qual = (char*) malloc(allocatedReadLen*sizeof(char));

		c = (char) getc(readsFile);
		while (c != '\n' && !feof(readsFile)) {
			if (readLen >= allocatedReadLen) {
				//printf("Loading reads: readLen = %u, allocReadLen = %u, reallocing...", readLen, allocatedReadLen);
				allocatedReadLen <<= 1;
				read->seq = (char*) realloc(read->seq, allocatedReadLen*sizeof(char));
				read->rc = (char*) realloc(read->rc, allocatedReadLen*sizeof(char));
				//read->qual = (char*) realloc(read->qual, allocatedReadLen*sizeof(char));
				//printf(" Success. \n");
			}
			read->seq[readLen] = nt4_table[(int) c];
			c = (char) getc(readsFile);
			readLen++;
		}
		read->len = readLen;
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
		int qualLen = 0;
		c = (char) getc(readsFile);
		while(c != '\n' && !feof(readsFile)) {
			if(qualLen <= readLen) {
				//read->qual[qualLen] = c;
			}
			qualLen++;
			c = (char) getc(readsFile);
		}
		if(qualLen != readLen) {
			printf("Error: The number of quality score symbols does not match the length of the read sequence.\n");
			exit(1);
		}
		//read->qual[qualLen] = '\0';

		// compute the reverse complement
		for(int i = 0; i < read->len; i++) {
			read->rc[read->len-1-i] = nt4_complement[(int)read->seq[i]];
		}
		reads->count++;
		if(reads->count > MAX_NUM_READS) {
			printf("Warning: The number of reads in %s exceeds the maximum number of reads allowed (= %d). "
					"The remaining read sequences will be ignored.\n", readsFname, MAX_NUM_READS);
			break;
		}
	}
	fclose(readsFile);
	return reads;
}

char *strdup(const char *str) {
    int n = strlen(str) + 1;
    char *dup = malloc(n * sizeof(char));
    if(dup)
    {
        strcpy(dup, str);
    }
    return dup;
}

// assumes that reads were generated with wgsim
void parse_read_mapping(char* read_name, unsigned int* ref_pos_l, unsigned int* ref_pos_r, int* strand) {
    int token_index = 0;
    const char delimiters[] = "_";
    char* read_name_dup = strdup(read_name);
    char* token = strtok(read_name_dup, delimiters);
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
        token = strtok(NULL, delimiters);
        token_index++;
    }
    free(read_name_dup);
}

void free_reads(reads_t* reads) {
	for(int i = 0; i < reads->count; i++) {
		if(reads->reads[i].seq) free(reads->reads[i].seq);
		if(reads->reads[i].rc) free(reads->reads[i].rc);
		if(reads->reads[i].qual) free(reads->reads[i].qual);
	}
	free(reads->reads);
	free(reads);
}

void print_read(read_t* read) {
	printf("%s \n", read->name);
	//printf("FWD: ");
	for(int i = 0; i < read->len; i++) {
		printf("%c", iupacChar[(int) read->seq[i]]);
	} printf("\n");
	//printf("RC: ");
	//for(int i = 0; i < read->len; i++) {
		//printf("%c", iupacChar[(int) read->rc[i]]);
	//} printf("\n");
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