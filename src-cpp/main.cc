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
#include "mt64.h"

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

void compute_hash_diff_stats(ref_t& ref, const reads_t& reads, const index_params_t* params) {
	printf("**********Diff Stats**************\n");

	VectorU32 diff_hist(SIMHASH_BITLEN);
	VectorU32 minh_hist(params->h);
	uint32 n_unmapped = 0;

	uint32 n_checked = 0;

	#pragma omp parallel for
	for(uint32 i = 0; i < reads.reads.size(); i++) {
			read_t r = reads.reads[i];
			unsigned int pos_l, pos_r;
			int strand;
			parse_read_mapping(r.name.c_str(), &pos_l, &pos_r, &strand);
			seq_t true_pos = pos_l - 1;

			MapPos2Window::iterator v;
			if((v = ref.windows_by_pos.find(true_pos)) == ref.windows_by_pos.end()) {
							n_unmapped++;
							continue;
			}
			ref_win_t ref_window = v->second;

			if(params->alg == MINH) {
					// count the number of minh values shared between the read and its window
					uint32 n_minh_shared = 0;
					for(uint32 h = 0; h < params->h; h++) {
						minhash_t minh = r.minhashes[h];
						for(uint32 h_ref = 0; h_ref < params->h; h_ref++) {
								if(ref_window.minhashes[h_ref] == minh) {
										n_minh_shared++;
										break;
								}
						}
					}

					#pragma omp critical
					{
						if(n_minh_shared == 0 && n_checked < 5) {
								printf("%s \n", r.name.c_str());
								print_read(&r);
								for(uint32 p = 0; p < params->ref_window_size; p++) {
										printf("%c", iupacChar[(int) ref.seq[true_pos + p]]);
								} printf("\n");
								printf("%d %llu \n", strand, true_pos);

								for(uint32 h = 0; h < params->h; h++) {
										printf("%u %u \n", r.minhashes[h], ref_window.minhashes[h]);
								}
								n_checked++;
						}
					}

					#pragma omp atomic
					minh_hist[n_minh_shared]++;

			} else {
					int d = hamming_dist(ref_window.simhash, r.simhash);
					#pragma omp atomic
					diff_hist[d]++;
			}
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
		printf("Usage: ./srx [options] <simh|minh|sample> <ref.fa> <reads.fq> \n");
		exit(1);
	}

	index_params_t params;
	params.set_default_index_params();

	int compute_diff_stats = 0;
	int c;
	while ((c = getopt(argc-1, argv+1, "i:o:k:w:d:L:H:m:p:ONSt:h:")) >= 0) {
		switch (c) {
			case 'i': params.in_index_fname = std::string(optarg); break;
			case 'o': params.out_index_fname = std::string(optarg); break;
			case 'k': params.k = atoi(optarg); break;
			case 'm': params.m = atoi(optarg); break;
			case 'p': params.p = atoi(optarg); break;
			case 'h': params.h = atoi(optarg); break;
			case 'w': params.ref_window_size = atoi(optarg); break;
			case 'd': params.max_hammd = atoi(optarg); break;
			case 'L': params.min_freq = atof(optarg); break;
			case 'H': params.max_freq = atof(optarg); break;
			case 't': params.n_threads = atoi(optarg); break;
			case 'O': params.kmer_type = OVERLAP; break;
			case 'N': params.kmer_type = NON_OVERLAP; break;
			case 'S': compute_diff_stats = 1; break;
			default: return 0;
		}
	}

	printf("**********SRX**************\n");
	printf("SRX Parameters: \n");
	printf("k = %u \n", params.k);
	printf("m = %u \n", params.m);
	printf("p = %u \n", params.p);
	printf("max_hammd = %u \n", params.max_hammd);


	if (params.kmer_type == NON_OVERLAP) {
		params.kmer_dist = params.k;
	}

	// generate the sparse k-mer indices
	if(params.kmer_type == SPARSE && params.in_index_fname.empty()) {
		VectorU32 idx(params.ref_window_size);
		params.sparse_kmers.resize(params.m * params.k);

		for(uint32 i = 0; i < params.m; i++) {
			for(uint32 k = 0; k < params.ref_window_size; k++) {
				idx[k] = k;
			}
			// pick random k *close-by* indices from the read
			const int32_t offset = i*params.k;
			const int32_t start = rand_range(params.ref_window_size - params.max_range);
			params.sparse_kmers[offset] = start;
			uint32 cnt = 1;
			uint32 len = params.max_range;
			while(cnt < params.k) {
				int j = rand_range(len) + 1; // exclude 0
				params.sparse_kmers[offset + cnt] = start + idx[j];
				idx[j] = idx[len-1];
				cnt++;
				len--;
			}
		}
	}

	// prepare the input and output files
	std::string histFname;
	std::string sparseFname;
	std::string hashFname;
	if(!params.in_index_fname.empty() || !params.out_index_fname.empty()) {
		std::string fname = params.in_index_fname;
		if(params.in_index_fname.empty()) {
			fname = params.out_index_fname;
		}
		histFname += fname + std::string(".hst");
		sparseFname += fname + std::string(".sparse");
		hashFname += fname + std::string(".hash");
	}

	if (strcmp(argv[1], "sample") == 0) {
		printf("LSH Algorithm: sampling \n");
		params.alg = SAMPLE;

		// 1. load the reference
		ref_t ref;
		fasta2ref(argv[optind+1], ref);

		// 2. compute valid reference windows
		generate_ref_windows(ref, &params);
		printf("Total number of valid windows: %zu\n", ref.windows_by_pos.size());

		// 3. load the reads
		reads_t reads;
		fastq2reads(argv[optind+2], reads);

		// 4. map by sampling
		align_reads_sampling(ref, reads, &params);

	} else if (strcmp(argv[1], "simh") == 0) { // SIMHASH
		printf("LSH Algorithm: simhash \n");
		params.alg = SIMH;

		// 1. index the reference
		ref_t ref;
		if(!params.in_index_fname.empty()) { // load the index
			load_ref_idx(params.in_index_fname.c_str(), ref);
			load_perm(sparseFname.c_str(), params.sparse_kmers);
		} else {
			index_ref_lsh(argv[optind+1], &params, ref);
		}
		// store the index
		if(!params.out_index_fname.empty()) {
			store_ref_idx(params.out_index_fname.c_str(), ref);
			if (params.kmer_type == SPARSE) {
				store_perm(sparseFname.c_str(), params.sparse_kmers);
			}
		}

		// 2. index the reads
		reads_t reads;
		index_reads_lsh(argv[optind+2], ref, &params, reads);

		if(compute_diff_stats) {
			compute_hash_diff_stats(ref, reads, &params);
		}

		// 3. map the hashes
		//align_reads_lsh(ref, reads, &params);

	} else if (strcmp(argv[1], "minh") == 0) { // MINHASH
		printf("LSH Algorithm: minhash \n");
		params.alg = MINH;

		// generate the independent hash functions
		params.rand_hash_pads.resize(params.h);
		for(uint32 i = 0; i < params.h; i++) {
			params.rand_hash_pads[i] = genrand64_int64();
			//printf("Hash %llx \n", params.rand_hash_pads[i]);
		}

		// 1. index the reference
		ref_t ref;
		if(!params.in_index_fname.empty()) {
			load_ref_idx(params.in_index_fname.c_str(), ref);
			load_perm(sparseFname.c_str(), params.sparse_kmers);
		} else {
			index_ref_lsh(argv[optind+1], &params, ref);
		}
		// store the index
		if(!params.out_index_fname.empty()) {
			store_ref_idx(params.in_index_fname.c_str(), ref);
			if (params.kmer_type == SPARSE) {
				store_perm(sparseFname.c_str(), params.sparse_kmers);
			}
			// TODO: store the hash functions
		}

		// 2. index the reads
		reads_t reads;
		index_reads_lsh(argv[optind+2], ref, &params, reads);

		if(compute_diff_stats) {
			compute_hash_diff_stats(ref, reads, &params);
		}

		// 3. map the hashes
		align_reads_minhash(ref, reads, &params);
	} else {
		printf("Usage: ./srx [options] <simh|minh|sample> <ref.fa> <reads.fq> \n");
		exit(1);
	}

	return 0;
}
