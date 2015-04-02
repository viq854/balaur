#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>

#include "types.h"
#include "index.h"
#include "align.h"
#include "io.h"


void compute_hash_diff_stats(ref_t& ref, const reads_t& reads, const index_params_t* params) {
//	printf("**********Diff Stats**************\n");
//
//	VectorU32 diff_hist(SIMHASH_BITLEN);
//	VectorU32 minh_hist(params->h);
//	uint32 n_unmapped = 0;
//
//	uint32 n_checked = 0;
//
//	#pragma omp parallel for
//	for(uint32 i = 0; i < reads.reads.size(); i++) {
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
//	}
//
//	if(params->alg == MINH) {
//		for(uint32 i = 0; i < params->h; i++) {
//			printf("%d: %d\n", i, minh_hist[i]);
//		}
//		printf("unmapped: %d\n", n_unmapped);
//	}
//	else {
//		for(uint32 i = 0; i < SIMHASH_BITLEN; i++) {
//			printf("%d: %d\n", i, diff_hist[i]);
//		}
//	}
}

int main(int argc, char *argv[]) {
	if (argc < 4) {
		printf("Usage: ./rx [options] <index|align|cluster> <ref.fa> <reads.fq> \n");
		exit(1);
	}

	index_params_t params;
	params.set_default_index_params();

	int c;
	while ((c = getopt(argc-1, argv+1, "i:o:w:k:h:L:H:T:b:p:l:g:St:I:m:s:d:v:N:")) >= 0) {
		switch (c) {
			case 'i': params.in_index_fname = std::string(optarg); break;
			case 'o': params.out_index_fname = std::string(optarg); break;
			case 'w': params.ref_window_size = atoi(optarg); break;
			case 'k': params.k = atoi(optarg); break;
			case 'v': params.k2 = atoi(optarg); break;
			case 'h': params.h = atoi(optarg); break;
			case 'H': params.max_count = atoi(optarg); break;
			case 'L': params.min_count = atoi(optarg); break;
			case 'T': params.n_tables = atoi(optarg); break;
			case 'b': params.sketch_proj_len = atoi(optarg); break;
			case 'p': params.n_buckets_pow2 = atoi(optarg); break;
			case 'l': params.bucket_entry_coverage = atoi(optarg); break;
			case 'g': params.contig_gap = atoi(optarg); break;
			case 't': params.n_threads = atoi(optarg); break;
			case 'I': params.hit_collection_interval = atoi(optarg); break;
			case 'm': params.min_n_hits = atoi(optarg); break;
			case 's': params.bucket_size = atoi(optarg); break;
			case 'd': params.dist_best_hit = atoi(optarg); break;
			case 'N': params.n_top_buckets_search = atoi(optarg); break;
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
	if (strcmp(argv[1], "index") == 0) {
		printf("Mode: Indexing \n");
		params.alg = MINH; // only min-hash enabled for now
		params.set_kmer_hash_function();
		params.set_minhash_hash_function();
		params.set_minhash_sketch_hash_function();
		params.generate_sparse_sketch_projections();

		// index the reference
		ref_t ref;
		index_ref_lsh(argv[optind+1], &params, ref);

		// store the index
		store_index_ref_lsh(argv[optind+1], &params, ref);

	} else if (strcmp(argv[1], "align") == 0) {
		printf("Mode: Alignment \n");
		params.alg = MINH; // only min-hash enabled for now
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

		// 3. map the hashes
		align_reads_minhash(ref, reads, &params);

	} else if (strcmp(argv[1], "cluster") == 0) {

	} else {
		printf("Usage: ./rx [options] <align|cluster> <ref.fa> <reads.fq> \n");
		exit(1);
	}

	return 0;
}
