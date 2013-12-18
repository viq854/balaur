#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "io.h"

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
		int allocatedReadLen = READ_LENGTH;
		read->seq = (char*) malloc(allocatedReadLen*sizeof(char));
		read->rc = (char*) malloc(allocatedReadLen*sizeof(char));
		//read->qual = (char*) malloc(allocatedReadLen*sizeof(char));

		c = (char) getc(readsFile);
		while (c != '\n' && !feof(readsFile)) {
			if (readLen >= allocatedReadLen) {
				//printf("Loading reads: readLen = %u, allocReadLen = %u, reallocing...", readLen, allocatedReadLen);
				allocatedReadLen += READ_LENGTH;
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
	printf("READ %s \n", read->name);
	printf("FWD: ");
	for(int i = 0; i < read->len; i++) {
		printf("%c", iupacChar[(int) read->seq[i]]);
	} printf("\n");
	printf("RC: ");
	for(int i = 0; i < read->len; i++) {
		printf("%c", iupacChar[(int) read->rc[i]]);
	} printf("\n");
}