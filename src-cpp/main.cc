/*
 * Program TOTORO
 * Victoria Popic (viq@stanford.edu) 2013-2015
 *
 * MIT License
 *
 * Copyright (c) 2015 Victoria Popic.
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include "index.h"
#include "align.h"

void print_usage(index_params_t* params) {
	printf("Usage: ./totoro [options] <index|align> <ref.fa> <reads.fq> \n");
	printf("Hashing options:\n\n");
	printf("       -h        number of hash functions for MinHash fingerprint construction (i.e. fingerprint length) [%d]\n", params->h);
	printf("       -T        number of hash tables [%d]\n", params->n_tables);
	printf("       -k        length of the sequence kmers [%d]\n", params->k);
	printf("       -b        length of the fingerprint projections [%d]\n", params->sketch_proj_len);
	printf("\nIndex-only options:\n\n");
	printf("       -w        length of the reference windows to hash (should be set to the expected read length for optimal results) [%d]\n", params->ref_window_size);
	printf("       -H        upper bound on kmer occurrence in the reference [%llu]\n", params->max_count);
	printf("       -s        initially allocated hash table bucket size [%d]\n", params->bucket_size);
	printf("\nAlignment-only options:\n\n");
	printf("       -m        minimum required number of buckets shared between a reference window and the read for a contig to be examined [%d]\n", params->min_n_hits);
	printf("       -N        maximum distance from the best number of shared buckets found for a contig to be examined [%d]\n", params->dist_best_hit);
	printf("       -v        length k2 of kmers counted during voting [%d]\n", params->k2);
	printf("       -P        disable the pre-computation of k2 hash values for the reference kmers (executed once per k2 length) [ON]\n");
	printf("       -n        number of initial inlier anchors to consider [%d]\n", params->n_init_anchors);
	printf("       -d        delta distance from the inlier anchor median considered close enough [%d]\n", params->delta_inlier);
	printf("       -x        delta multiplier for the second RANSAC pass (i.e. number of deltas away the second position must be from the first pass median) [%d]\n", params->delta_x);
	printf("       -c        cutoff minimum number of inlier votes [dynamic]\n");
	printf("       -S        disable votes scaling [ON]\n");
	printf("\nOther options:\n\n");
	printf("       -t        number of threads [%d]\n", params->n_threads);
}

int main(int argc, char *argv[]) {
	index_params_t params;
	params.set_default_index_params();

	if (argc < 4) {
		print_usage(&params);
		exit(1);
	}
	int c;
	while ((c = getopt(argc-1, argv+1, "i:o:w:k:h:L:H:T:b:p:l:t:m:s:d:v:PN:n:c:Sx:")) >= 0) {
		switch (c) {
			case 'h': params.h = atoi(optarg); break;
			case 'T': params.n_tables = atoi(optarg); break;
			case 'k': params.k = atoi(optarg); break;
			case 'b': params.sketch_proj_len = atoi(optarg); break;
			case 'w': params.ref_window_size = atoi(optarg); break;
			case 'p': params.n_buckets_pow2 = atoi(optarg); break;
			case 's': params.bucket_size = atoi(optarg); break;
			case 'l': params.bucket_entry_coverage = atoi(optarg); break;
			case 'H': params.max_count = atoi(optarg); break;
			case 'L': params.min_count = atoi(optarg); break;
			case 'm': params.min_n_hits = atoi(optarg); break;
			case 'N': params.dist_best_hit = atoi(optarg); break;
			case 'v': params.k2 = atoi(optarg); break;
			case 'P': params.precomp_k2 = false; break;
			case 'S': params.enable_scale = false; break;
			case 'n': params.n_init_anchors = atoi(optarg); break;
			case 'd': params.delta_inlier = atoi(optarg); break;
			case 'x': params.delta_x = atoi(optarg); break;
			case 'c': params.votes_cutoff = atoi(optarg); break;
			case 'i': params.in_index_fname = std::string(optarg); break;
			case 'o': params.out_index_fname = std::string(optarg); break;
			case 't': params.n_threads = atoi(optarg); break;
			default: return 0;
		}
	}

	printf("**********TOTORO**************\n");
	srand(1);
	params.n_buckets = pow(2, params.n_buckets_pow2);
	if (strcmp(argv[1], "index") == 0) {
		printf("Mode: Indexing \n");
		params.alg = MINH; // only minhash enabled for now
		params.set_kmer_hash_function();
		params.set_minhash_hash_function();
		params.set_minhash_sketch_hash_function();
		params.generate_sparse_sketch_projections();
		// index the reference and store
		ref_t ref;
		index_ref_lsh(argv[optind+1], &params, ref);
		store_index_ref_lsh(argv[optind+1], &params, ref);

	} else if (strcmp(argv[1], "align") == 0) {
		printf("Mode: Alignment \n");
		params.alg = MINH; // only minhash enabled for now
		params.set_kmer_hash_function();
		params.set_minhash_hash_function();
		params.set_minhash_sketch_hash_function();
		params.generate_sparse_sketch_projections();

		// 1. load the reference index
		ref_t ref;
		load_index_ref_lsh(argv[optind+1], &params, ref);

		// 2. index the reads
		reads_t reads;
		index_reads_lsh(argv[optind+2], ref, &params, reads);

		// 3. align
		align_reads_minhash(ref, reads, &params);
	} else if (strcmp(argv[1], "stats") == 0) {
		printf("Mode: STATS \n");

		// load the reference index
		ref_t ref;
		load_index_ref_lsh(argv[optind+1], &params, ref);
		store_ref_index_stats(argv[optind+1], ref, &params);
	} else {
		print_usage(&params);
		exit(1);
	}
	return 0;
}
