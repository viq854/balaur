/*
 * Program BALAUR
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
#include "sam.h"

void print_usage() {
	printf("Usage: ./balaur [options] <index|align> <ref.fa> <reads.fq> \n");
	printf("Hashing options:\n\n");
	printf("       -h        number of hash functions for MinHash fingerprint construction (i.e. fingerprint length) [%d]\n", params->h);
	printf("       -T        number of hash tables [%d]\n", params->n_tables);
	printf("       -k        length of the sequence kmers [%d]\n", params->k);
	printf("       -b        length of the fingerprint projections [%d]\n", params->sketch_proj_len);
	printf("\nIndex-only options:\n\n");
	printf("       -w       length of the reference windows to hash (should be set to the expected read length for optimal results) [%d]\n", params->ref_window_size);
	printf("       -H       upper bound on kmer occurrence in the reference [%llu]\n", params->max_count);
	printf("       -s       initially allocated hash table bucket size [%d]\n", params->bucket_size);
	printf("\nAlignment-only options:\n\n");
	printf("       -v        length k2 of kmers counted during voting [%d]\n", params->k2);
	printf("       -m      	 candidate contig filtering: min required number of index buckets shared with the read [%d]\n", params->min_n_hits);
	printf("       -N        candidate contig filtering: max distance from the best found number of buckets shared with a contig [%d]\n", params->dist_best_hit);
	printf("       -L        load precomputed candidate contigs [%d]\n", params->k2);
	printf("       -z        precomputed candidate contigs file (store/load) [%d]\n", params->k2);
	printf("       -d        votes array convolution radius  [%d]\n", params->delta_inlier);
	printf("       -x        min distance between separate mapping hits  [%d]\n", params->delta_x);
	printf("       -c        cutoff on the min number of votes [auto]\n");
	printf("       -f        mapq scaling factor  [auto]\n");
	printf("       -I        voting kmer sampling factor [%d]\n", params->sampling_intv);
	printf("       -Z        read batch size (number of reads loaded/processed at a time) [%d]\n", READ_BATCH_SIZE);
	printf("\nPrivacy-related options:\n\n");
	printf("       -V       enable vanilla mode (non-cryptographic hashing, no repeat filtering) \n");
	printf("       -B       voting kmer discretized position range  [%d]\n", params->bin_size);
	printf("       -S       voting task size: number of contigs per read encrypted with same keys [%d]\n", params->batch_size);
	printf("       -M       enable masking kmers neighboring repeats (default: only repeats are masked) \n");
	printf("\nOther options:\n\n");
	printf("       -t        number of threads, index/voting computation [%d]\n", params->n_threads);
}

index_params_t* params;
int main(int argc, char *argv[]) {
	params = new index_params_t();
	params->set_default_index_params();
	uint32 read_batch_size = READ_BATCH_SIZE;

	if (argc < 4) {
		print_usage();
		exit(1);
	}
	int c;
	while ((c = getopt(argc-1, argv+1, "t:w:k:h:H:T:b:p:m:s:d:v:N:c:x:Lf:z:I:Z:S:B:MV")) >= 0) {
		switch (c) {
			case 't': params->n_threads = atoi(optarg); break;
			case 'h': params->h = atoi(optarg); break;
			case 'T': params->n_tables = atoi(optarg); break;
			case 'k': params->k = atoi(optarg); break;
			case 'b': params->sketch_proj_len = atoi(optarg); break;
			case 'w': params->ref_window_size = atoi(optarg); break;
			case 'p': params->n_buckets_pow2 = atoi(optarg); break;
			case 's': params->bucket_size = atoi(optarg); break;
			case 'H': params->max_count = atoi(optarg); break;
			case 'm': params->min_n_hits = atoi(optarg); break;
			case 'N': params->dist_best_hit = atoi(optarg); break;
			case 'L': params->load_mhi = false; break;
			case 'z': params->precomp_contig_file_name = std::string(optarg); break;
			case 'v': params->k2 = atoi(optarg); break;
			case 'd': params->delta_inlier = atoi(optarg); break;
			case 'x': params->delta_x = atoi(optarg); break;
			case 'c': params->votes_cutoff = atoi(optarg); break;
			case 'f': params->mapq_scale_x = atoi(optarg); break;
			case 'I': params->sampling_intv = atoi(optarg); break;
			case 'V': params->vanilla = true; break;
			case 'S': params->batch_size = atoi(optarg); break;
			case 'B': params->bin_size = atoi(optarg); break;
			case 'M': params->mask_repeat_nbrs = true; break;
			case 'Z': read_batch_size = atoi(optarg); break;
			default: return 0;
		}
	}
	srand(1);
	params->set_kmer_hash_function();
	params->set_minhash_hash_function();
	params->set_minhash_sketch_hash_function();
	params->generate_sparse_sketch_projections();
	params->set_bin_shuffle();
	if(params->vanilla) {
		printf("Note: VANILLA mode activated (privacy-related parameter settings will be ignored) \n");
		params->kmer_hashing_alg = CITY_HASH64;
		params->batch_size = INT_MAX;
		params->bin_size = 1;
		params->mask_repeat_nbrs = false;
		if(params->ref_window_size > 350) { // todo: expose params
			printf("Switching to *monolith* mode for long reads.\n");
			params->monolith = true;
		}
	}
		
	printf("**********BALAUR**************\n");
	params->n_buckets = pow(2, params->n_buckets_pow2);
	if (strcmp(argv[1], "index") == 0) {
		ref_t ref;
		index_ref_lsh(argv[optind+1], params, ref);
		store_index_ref_lsh(argv[optind+1], params, ref);
	} else if (strcmp(argv[1], "align") == 0) {
		// load the reference, index, and auxiliary ref structures
		ref_t ref;
		load_index_ref_lsh(argv[optind+1], params, ref);
		if(!load_kmer2_hashes(argv[optind+1], ref, params)) {
			printf("Precomputing voting kmer reference hashes; this only needs to run once for this -v value and encryption algorithm.\n");
			compute_store_kmer2_hashes(argv[optind+1], ref, params);
		}
		if(!params->vanilla) {
			if(!load_repeat_info(argv[optind+1], ref, params)) {
				printf("Precomputing kmer reference repeats; this only needs to run once for this -v value.\n");
				compute_store_repeat_info(argv[optind+1], ref, params);
			}
 		}

		// load/store contig candidates (optional)	
		precomp_contig_io_t contig_io;
		if(params->load_mhi) {
			if(params->precomp_contig_file_name.size() != 0) {
				contig_io.open_file(params->precomp_contig_file_name.c_str(), STORE);
			}
		} else {
			contig_io.open_file(params->precomp_contig_file_name.c_str(), LOAD);
		}
		
		// will store alignment results
		sam_writer_t sam_io;
		sam_io.open_file(argv[optind+2]);
	
		// load all reads as a single batch
                // reads_t reads;
                // fastq2reads(argv[optind+2], reads);

		// read batch processing
                fastq_reader_t reader;
                reader.open_file(std::string(argv[optind+2]));	
		while(true) {
			printf(">>>>>>>Loading read batch.....\n");
			reads_t reads;
			if(!reader.load_next_read_batch(reads, read_batch_size)) break;
			//ref_t ref;
                        //load_index_ref_lsh(argv[optind+1], params, ref);
			balaur_main(argv[optind+1], ref, reads, contig_io, sam_io);
			reads.free();
		}
		printf(">>>>>>>Alignment done!\n");
		reader.close_file();
		contig_io.close_file();
		sam_io.close_file();
	} /*else if (strcmp(argv[1], "stats") == 0) {
		printf("Mode: STATS \n");
		//ref_t ref;
		//load_ref_idx(argv[optind+1], ref, &params);
		//store_ref_idx_flat(argv[optind+1], ref, &params);

		ref_t ref;
		fasta2ref(argv[optind+1], ref);
		load_kmer2_hashes(argv[optind+1], ref, params);
		compute_store_repeat_info(argv[optind+1], ref, params);
		//compute_store_repeat_local(argv[optind+1], ref, &params);		
		//load_repeat_info(argv[optind+1], ref, params);
		//bin_repeat_stats(argv[optind+1], params, ref);
		//compute_store_kmer2_hashes(argv[optind+1], ref, params);
		//ref_kmer_repeat_stats(argv[optind+1], &params, ref);
		//ref_kmer_fingerprint_stats(argv[optind+1], &params, ref);
		//load_index_ref_lsh(argv[optind+1], &params, ref);
		//store_ref_index_stats(argv[optind+1], ref, &params);
		//kmer_stats(argv[optind+1]);
	} */ else {
		std::cout<< "No valid command specified.\n";
		print_usage();
		exit(1);
	}
	return 0;
}
