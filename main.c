#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include "index.h"
#include "align.h"
#include "hash.h"
#include "mt64.h"

void set_default_index_params(index_params_t* params) {
	params->p = 1;
	params->k = 64;
	params->m = 2;
	params->h = 64;
	params->min_freq = 0.000001;
	params->max_freq = 0.6;
	params->ref_window_size = 100;
	params->max_hammd = 20;
}

void generate_reads(char* fname) {
	FILE* readsFile = (FILE*) fopen(fname, "w");
	if (readsFile == NULL) {
		printf("Cannot open reads file!\n");
		exit(1);
	}
	
	int readlen = 100;
	int readnum = 100;
	int pos = 0;
	for(int i = 0; i < readnum; i++) {
		fprintf(readsFile, "@r%d\n", i);
		for(int j = 0; j < readlen; j++) {
			if(j == pos) {
				fprintf(readsFile, "%c", 'C');
			} else {
				fprintf(readsFile, "%c", 'A');
			}
		}
		fprintf(readsFile, "\n+\n");
		for(int j = 0; j < readlen; j++) {
			fprintf(readsFile, "%d", 2);		
		}
		fprintf(readsFile, "\n");
		pos++;
		if(pos == readnum) {
			pos = 0;
		}
	}
	fclose(readsFile);
}


int rand_range(int n) {
	int r, rand_max = RAND_MAX - (RAND_MAX % n);
	while ((r = rand()) >= rand_max);
	return r / (rand_max / n);
}


int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: srx [options] <ref.fa> <reads.fq> \n");
		exit(1);
	}
	index_params_t* params = (index_params_t*) calloc(1, sizeof(index_params_t));
	set_default_index_params(params);

	int c;
	while ((c = getopt(argc-1, argv+1, "I:k:w:d:L:H:")) >= 0) {
		switch (c) {
			
			case 'k': params->k = atoi(optarg); break;
			case 'w': params->ref_window_size = atoi(optarg); break;
			case 'd': params->max_hammd = atoi(optarg); break;
			case 'L': params->min_freq = atof(optarg); break;
			case 'H': params->max_freq = atof(optarg); break;
			default: return 0;
		}
	}
	params->hist_size = KMER_HIST_SIZE16;
	if(params->k > CHARS_PER_SHORT) {
		params->hist_size = KMER_HIST_SIZE32;
	}
		
	// generate the sparse k-mer indices
	int idx[100] = { 0 };
	for(int i = 0; i < 100; i++) {
		idx[i] = i;
	}
	params->sparse_kmers = (int*) malloc(params->m*params->k*sizeof(int));
	for(int i = 0; i < params->m; i++) {
		for(int k = 0; k < 100; k++) {
			idx[k] = k;
        }
		int offset = i*params->k;
		// pick random k indices from the read
		int cnt = 0;
		int len = params->ref_window_size; 
		while(cnt < params->k) {		
			int j = rand_range(len);						
			params->sparse_kmers[offset+cnt] = idx[j];
			idx[j] = idx[len-1];
			cnt++;
			len--;			
		}
	}
	
	params->rand_hash_pads = (simhash_t*) malloc(params->h*sizeof(simhash_t));
	for(int i = 0; i < params->h; i++) {
		params->rand_hash_pads[i] = genrand64_int64();
	}

	// 1. ref
	ref_t* ref;
	//index_ref_windows(argv[1], params, &ref);
	index_ref_simhash(argv[1], params, &ref);
	
	// save index to file
	char* idxFname  = (char*) malloc(strlen(argv[1]) + 12);
	char* histFname  = (char*) malloc(strlen(argv[1]) + 12);
	char* permFname  = (char*) malloc(strlen(argv[1]) + 12);
	char* padFname  = (char*) malloc(strlen(argv[1]) + 12);
	sprintf(idxFname, "%s.idx_sp%d.%d", argv[1], params->k, params->m);
	sprintf(histFname, "%s.hst_sp%d.%d", argv[1], params->k, params->m);
	sprintf(permFname, "%s.perm%d.%d", argv[1], params->k, params->m);
	sprintf(padFname, "%s.pad%d", argv[1], params->h);
	//store_ref_idx(ref, idxFname);	
	//ref = load_ref_idx(idxFname);
	//store_perm(params->sparse_kmers, params->m*params->k, permFname);
	//params->sparse_kmers = load_perm(params->m*params->k, permFname);
	
	params->max_count = (uint64_t) (params->max_freq*ref->num_windows);
	
	free(histFname);
	free(idxFname);
	free(permFname);
	free(padFname);
	
	// 2. reads
	reads_t* reads;
	index_reads_simhash(argv[2], ref, params, &reads);
	//reads = fastq2reads(argv[2]);
	
	// 3. map
	//align_reads_sampling(ref, reads, params);
	align_reads_minhash(ref, reads, params);

	free(params);
	return 0;
}
