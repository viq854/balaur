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
	params->k = 16;
	params->hist_size = KMER_HIST_SIZE16;
	params->m = 10;
	params->max_range = params->k + 10;
	params->p = 1;
	params->h = 64;
	params->min_freq = 0.000001;
	params->max_freq = 0.6;
	params->ref_window_size = 100;
	params->max_hammd = 10;
	params->kmer_type = SPARSE;
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
	if (argc < 4) {
		printf("Usage: ./srx [options] <simh|minh|sample> <ref.fa> <reads.fq> \n");
		exit(1);
	}

	index_params_t* params = (index_params_t*) calloc(1, sizeof(index_params_t));
	set_default_index_params(params);

	int c;
	while ((c = getopt(argc-1, argv+1, "i:o:k:w:d:L:H:m:p:ON")) >= 0) {
		switch (c) {
			case 'i': params->in_idx_fname = optarg; break;
			case 'o': params->out_idx_fname = optarg; break;
			case 'k': params->k = atoi(optarg); break;
			case 'm': params->m = atoi(optarg); break;
			case 'p': params->p = atoi(optarg); break;
			case 'w': params->ref_window_size = atoi(optarg); break;
			case 'd': params->max_hammd = atoi(optarg); break;
			case 'L': params->min_freq = atof(optarg); break;
			case 'H': params->max_freq = atof(optarg); break;
			case 'O': params->kmer_type = OVERLAP; break;
			case 'N': params->kmer_type = NON_OVERLAP; break;
			default: return 0;
		}
	}

	// select the size of the reference kmer histogram table needed
	if(params->k > CHARS_PER_SHORT) {
		params->hist_size = KMER_HIST_SIZE32;
	}

	// generate the sparse k-mer indices
	if(params->in_idx_fname == NULL && params->kmer_type == SPARSE) {
		int* idx = (int*) malloc(params->ref_window_size*sizeof(int));
		params->sparse_kmers = (int*) malloc(params->m*params->k*sizeof(int));
		for(int i = 0; i < params->m; i++) {
			for(int k = 0; k < params->ref_window_size; k++) {
				idx[k] = k;
			}
			// pick random k *close-by* indices from the read
			int offset = i*params->k;
			int start = rand_range(params->ref_window_size - params->max_range);
			params->sparse_kmers[offset] = start;
			int cnt = 1;
			int len = params->max_range;
			while(cnt < params->k) {
				int j = rand_range(len) + 1; // exclude 0
				params->sparse_kmers[offset + cnt] = start + idx[j];
				idx[j] = idx[len-1];
				cnt++;
				len--;
			}
		}
		free(idx);

#if (DEBUG)
		for(int i = 0; i < params->m; i++) {
			printf("m = %d \n", i);
			for(int j = 0; j < params->k; j++) {
				printf(" %d ", params->sparse_kmers[i*params->k + j]);
			}
			printf("\n");
		}
#endif
	}

	// prepare the input and output files
	char* histFname;
	char* sparseFname;
	char* hashFname;
	if(params->in_idx_fname || params->out_idx_fname) {
		const char* fname = params->in_idx_fname;
		if(params->in_idx_fname) {
			fname = params->out_idx_fname;
		}
		histFname  = (char*) malloc(strlen(fname) + 12);
		sparseFname  = (char*) malloc(strlen(fname) + 12);
		hashFname  = (char*) malloc(strlen(fname) + 12);
		sprintf(histFname, "%s.hst", fname);
		sprintf(sparseFname, "%s.sparse", fname);
		sprintf(hashFname, "%s.hash", fname);
	}

	if (strcmp(argv[1], "simh") == 0) { // SIMHASH
		printf("LSH Algorithm: simhash \n");

		params->alg = SIMH;

		// 1. index the reference
		ref_t* ref;
		// load the index
		if(params->in_idx_fname) {
			ref = load_ref_idx(params->in_idx_fname);
			params->sparse_kmers = load_perm(params->m*params->k, sparseFname);
		} else {
			index_ref_lsh(argv[optind+1], params, &ref);
		}
		// store the index
		if(params->out_idx_fname) {
			store_ref_idx(ref, params->in_idx_fname);
			if (params->kmer_type == SPARSE) {
				store_perm(params->sparse_kmers, params->m*params->k, sparseFname);
			}
		}

		// 2. index the reads
		reads_t* reads;
		index_reads_lsh(argv[optind+2], ref, params, &reads);

		// 3. map the hashes
		align_reads_simhash(ref, reads, params);

	} else if (strcmp(argv[1], "minh") == 0) { // MINHASH
		printf("LSH Algorithm: minhash \n");

		params->alg = MINH;

		// 1. generate the independent hash functions
		params->rand_hash_pads = (hash_t*) malloc(params->h*sizeof(hash_t));
		for(int i = 0; i < params->h; i++) {
			params->rand_hash_pads[i] = genrand64_int64();
		}

		// 2. hash the reference
		ref_t* ref;
		// load the index
		if(params->in_idx_fname) {
			ref = load_ref_idx(params->in_idx_fname);
			params->sparse_kmers = load_perm(params->m*params->k, sparseFname);
		} else {
			index_ref_lsh(argv[optind+1], params, &ref);
		}
		// store the index
		if(params->out_idx_fname) {
			store_ref_idx(ref, params->in_idx_fname);
			store_perm(params->sparse_kmers, params->m*params->k, sparseFname);
			// TODO: store the hash functions
		}

		// 3. hash the reads
		reads_t* reads;
		index_reads_lsh(argv[optind+2], ref, params, &reads);

		// 4. map the hashes
		align_reads_minhash(ref, reads, params);

	} else if (strcmp(argv[1], "sample") == 0) {
		printf("LSH Algorithm: sampling \n");

		params->alg = SAMPLE;

		// 1. load the reference
		ref_t* ref = fasta2ref(argv[optind+1]);

		// 2. compute valid reference windows
		generate_ref_windows(ref, params);
		printf("Total number of valid windows: %llu\n", ref->num_windows);

		// 3. load the reads
		reads_t* reads = fastq2reads(argv[optind+2]);

		// 4. map by sampling
		align_reads_sampling(ref, reads, params);

	} else {
		printf("Usage: ./srx [options] <simh|minh|sample> <ref.fa> <reads.fq> \n");
		exit(1);
	}

	if(params->in_idx_fname || params->out_idx_fname) {
		free(histFname);
		free(sparseFname);
		free(hashFname);
	}
	free(params);
	return 0;
}
