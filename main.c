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

void set_default_index_params(index_params_t* params) {
	params->k = 4;
	params->min_freq = 0.001;
	params->max_freq = 0.6;
	params->ref_window_size = 6;
	params->max_hammd = 7;
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

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: srx [options] <ref.fa> <reads.fq> \n");
		exit(1);
	}
	index_params_t* params = (index_params_t*) calloc(1, sizeof(index_params_t));
	set_default_index_params(params);

	int c;
	while ((c = getopt(argc-1, argv+1, "k:w:d:")) >= 0) {
		switch (c) {
			case 'k': params->k = atoi(optarg); break;
			case 'w': params->ref_window_size = atoi(optarg); break;
			case 'd': params->max_hammd = atoi(optarg); break;
			case 'minf': params->min_freq = atof(optarg); break;
			case 'maxf': params->max_freq = atof(optarg); break;
			default: return 0;
		}
	}
	
	//generate_reads(argv[1]);
	//index_reads(argv[1], params);
	
	ref_t* ref;
	index_ref(argv[1], params, &ref);
	reads_t* reads;
	index_reads(argv[2], ref, params, &reads);
	align_reads(ref, reads);
	
	/*for(uint64_t i = 0; i < 10; i++) {
		//uint64_t gray = (num >> 1) ^ num;
		//uint64_t num = i;
		//uint64_t mask;
	    //for (mask = num >> 1; mask != 0; mask = mask >> 1) {
	      //  num = num ^ mask;
	    //}
		uint64_t n = i; 
		uint64_t p = n;
		while (n >>= 1ULL) {
			p ^= n;
		}
		printf("%llu %llu \n",  i, p);
	}*/
	
	free(params);
	return 0;
}