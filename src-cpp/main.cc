#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <string.h>

#include "index.h"
#include "align.h"
#include "hash.h"


int rand_range(int n) {
	int r, rand_max = RAND_MAX - (RAND_MAX % n);
	while ((r = rand()) >= rand_max);
	return r / (rand_max / n);
}

void compute_hash_diff_stats(ref_t& ref, const reads_t& reads, const index_params_t* params) {
	printf("**********Diff Stats**************\n");

	VectorU32 diff_hist(SIMHASH_BITLEN);
	VectorU32 minh_hist(params->h);
	uint32 n_unmapped = 0;

	uint32 n_checked = 0;

	#pragma omp parallel for
	for(uint32 i = 0; i < reads.reads.size(); i++) {
//			read_t r = reads.reads[i];
//			unsigned int pos_l, pos_r;
//			int strand;
//			parse_read_mapping(r.name.c_str(), &pos_l, &pos_r, &strand);
//			seq_t true_pos = pos_l - 1;
//
//			MapPos2Window::iterator v;
//			if((v = ref.windows_by_pos.find(true_pos)) == ref.windows_by_pos.end()) {
//				print_read(&r);
//				n_unmapped++;
//				continue;
//			}
//			ref_win_t ref_window = v->second;
//
//			if(params->alg == MINH) {
//					// count the number of minh values shared between the read and its window
//					uint32 n_minh_shared = 0;
//					for(uint32 h = 0; h < params->h; h++) {
//						minhash_t minh = r.minhashes[h];
//						for(uint32 h_ref = 0; h_ref < params->h; h_ref++) {
//								if(ref_window.minhashes[h_ref] == minh) {
//										n_minh_shared++;
//										break;
//								}
//						}
//					}
//
//					#pragma omp critical
//					{
//						if(n_minh_shared == 0 && n_checked < 5) {
//								//printf("%s \n", r.name.c_str());
//								//print_read(&r);
//								//for(uint32 p = 0; p < params->ref_window_size; p++) {
//								///		printf("%c", iupacChar[(int) ref.seq[true_pos + p]]);
//								//} printf("\n");
//								//printf("%d %llu \n", strand, true_pos);
//
//								for(uint32 h = 0; h < params->h; h++) {
//									//printf("%u %u \n", r.minhashes[h], ref_window.minhashes[h]);
//								}
//								n_checked++;
//						}
//					}
//
//					#pragma omp atomic
//					minh_hist[n_minh_shared]++;
//
//			} else {
//					int d = hamming_dist(ref_window.simhash, r.simhash);
//					#pragma omp atomic
//					diff_hist[d]++;
//			}
	}

	if(params->alg == MINH) {
		for(uint32 i = 0; i < params->h; i++) {
			printf("%d: %d\n", i, minh_hist[i]);
		}
		printf("unmapped: %d\n", n_unmapped);
	}
	else {
		for(uint32 i = 0; i < SIMHASH_BITLEN; i++) {
			printf("%d: %d\n", i, diff_hist[i]);
		}
	}
}

int main(int argc, char *argv[]) {
	if (argc < 4) {
		printf("Usage: ./rx [options] <align|cluster> <ref.fa> <reads.fq> \n");
		exit(1);
	}

	index_params_t params;
	params.set_default_index_params();

	int compute_diff_stats = 0;
	int c;
	while ((c = getopt(argc-1, argv+1, "i:o:w:k:h:L:H:T:b:p:l:g:St:I:m:")) >= 0) {
		switch (c) {
			case 'i': params.in_index_fname = std::string(optarg); break;
			case 'o': params.out_index_fname = std::string(optarg); break;
			case 'w': params.ref_window_size = atoi(optarg); break;
			case 'k': params.k = atoi(optarg); break;
			case 'h': params.h = atoi(optarg); break;
			case 'H': params.max_count = atoi(optarg); break;
			case 'L': params.min_count = atoi(optarg); break;
			case 'T': params.n_tables = atoi(optarg); break;
			case 'b': params.sketch_proj_len = atoi(optarg); break;
			case 'p': params.n_buckets_pow2 = atoi(optarg); break;
			case 'l': params.bucket_entry_coverage = atoi(optarg); break;
			case 'g': params.contig_gap = atoi(optarg); break;
			case 't': params.n_threads = atoi(optarg); break;
			case 'S': compute_diff_stats = 1; break;
			case 'I': params.hit_collection_interval = atoi(optarg); break;
			case 'm': params.min_n_hits = atoi(optarg); break;
			default: return 0;
		}
	}

	printf("**********RX**************\n");
	printf("Configuration: \n");
	printf("w = %u \n", params.ref_window_size);
	printf("k = %u \n", params.k);
	printf("h = %u \n", params.h);
	printf("T = %u \n", params.n_tables);
	printf("b = %u \n", params.sketch_proj_len);
	printf("m = %u \n", params.min_n_hits);
	printf("n_buckets = %f \n", pow(2, params.n_buckets_pow2));


	srand(1);
	if (strcmp(argv[1], "align") == 0) {
		printf("Mode: Alignment \n");
		params.alg = MINH; // only min-hash enabled for now

		// prepare the input/output reference index files
		index_files_t index_files;
		index_files.prep_index_files(params.in_index_fname.empty() ? params.out_index_fname : params.in_index_fname);

		// set the initial kmer hash function (rolling hash)
		params.kmer_hasher = new CyclicHash(params.k, 32);

		// generate random hash functions for min-hash sketches
		for(uint32 h = 0; h < params.h; h++) {
			params.minhash_functions.push_back(rand_hash_function_t());
		}

		// generate random vector hash function sketch buckets
		params.sketch_proj_hash_func = rand_hash_function_t(params.n_buckets_pow2, params.sketch_proj_len);

		// generate the sparse sketch projections
		if(params.in_index_fname.empty()) {
			VectorU32 idx(params.h); // sketch length
			params.sketch_proj_indices.resize(params.sketch_proj_len * params.n_tables);

			for(uint32 i = 0; i < params.n_tables; i++) {
				for(uint32 k = 0; k < params.h; k++) {
					idx[k] = k;
				}
				// pick random indices from the sketch
				const int32_t offset = i*params.sketch_proj_len;
				const int32_t start = rand_range(params.h);
				params.sketch_proj_indices[offset] = start;
				uint32 cnt = 0;
				uint32 len = params.h;
				while(cnt < params.sketch_proj_len) {
					int j = rand_range(len); // exclude 0
					params.sketch_proj_indices[offset + cnt] = idx[j];
					idx[j] = idx[len-1];
					cnt++;
					len--;
				}
			}
		}

		// 1. index the reference
		ref_t ref;
		if(!params.in_index_fname.empty()) {
			load_ref_idx(params.in_index_fname.c_str(), ref);
		} else {
			index_ref_lsh(argv[optind+1], &params, ref);
		}
		// store the index
		if(!params.out_index_fname.empty()) {
			store_ref_idx(params.in_index_fname.c_str(), ref);
		}

		// 2. index the reads
		reads_t reads;
		index_reads_lsh(argv[optind+2], ref, &params, reads);

		if(compute_diff_stats) {
			compute_hash_diff_stats(ref, reads, &params);
		}

		// 3. map the hashes
		align_reads_minhash(ref, reads, &params);

	} else if (strcmp(argv[1], "cluster") == 0) {

	} else {
		printf("Usage: ./rx [options] <align|cluster> <ref.fa> <reads.fq> \n");
		exit(1);
	}

	return 0;
}
