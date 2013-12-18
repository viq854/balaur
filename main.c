#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <stdlib.h>
#include <getopt.h>
#include "index.h"


void set_default_index_params(index_params_t* params) {
	params->k = 2;
	params->s = 15;
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: srx [options] <reads.fq> \n");
		exit(1);
	}
	index_params_t* params = (index_params_t*) calloc(1, sizeof(index_params_t));
	set_default_index_params(params);

	int c;
	while ((c = getopt(argc-1, argv+1, "k:s:")) >= 0) {
		switch (c) {
			case 'k': params->k = atoi(optarg); break;
			case 's': params->s = atoi(optarg); break;
			default: return 0;
		}
	}
	index_reads(argv[1], params);
	free(params);
	return 0;
}