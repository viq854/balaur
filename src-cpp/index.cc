#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <limits.h>
#include "types.h"
#include "index.h"
#include "lsh.h"
#include "io.h"
#include "hash.h"
#include <fstream>

// --- Minhash ---

// index reference:
// - load fasta
// - generate valid windows (sliding window)
// - compute the fingerprint of each window
// - bucket into multiple hash tables using different projections
// - sort buckets

void ref_t::build_index(const char* fastaFname, const index_params_t* params) {
	double start_time = omp_get_wtime();
	fasta2ref(fastaFname, *this);
	load_freq_kmers(fastaFname, high_freq_kmer_bitmap, high_freq_kmer_trie, params->max_count);
	mark_freq_kmers(params);
	if(!load_valid_window_mask(fastaFname, *this, params)) {
		mark_windows_to_discard(params);
		store_valid_window_mask(fastaFname, *this, params);
	}
	printf("Total pre-processing time: %.2f sec\n", omp_get_wtime() - start_time);

	// initialize the hash tables
	mutable_index.per_table_buckets.resize(params->n_tables);
	for(uint32 t = 0; t < params->n_tables; t++) {
		buckets_t* buckets = &mutable_index.per_table_buckets[t];
		buckets->buckets_data_vectors.resize(params->n_buckets);
		buckets->per_thread_buckets_data_vectors.resize(params->n_threads);
		buckets->per_thread_bucket_sizes.resize(params->n_threads);
		for(uint32 tid = 0; tid < params->n_threads; tid++) {
			buckets->per_thread_buckets_data_vectors[tid].resize(params->n_buckets);
			buckets->per_thread_bucket_sizes[tid].resize(params->n_buckets, 0);
		}
	}

	// initialize additional per-thread storage
	std::vector<minhash_matrix_t> minhash_matrices(params->n_threads);
	std::vector<VectorMinHash> minhash_thread_vectors(params->n_threads);
	for(uint32 i = 0; i < params->n_threads; i++) {
		minhash_thread_vectors[i].resize(params->h);
		std::fill(minhash_thread_vectors[i].begin(), minhash_thread_vectors[i].end(), UINT_MAX);
	}

	// hash each valid window
	printf("Hashing reference windows... \n");
	uint32 n_valid_windows = 0;
	uint32 n_valid_hashes = 0;
	uint64 n_bucket_entries = 0;
	uint64 n_filtered = 0;

#if EXTERNAL_MEM_INDEX
#define MAX_NTABLES_NO_DISK 1024
	int file_nsync_points = 0;
	if(params->n_tables > MAX_NTABLES_NO_DISK) {
		file_nsync_points = 42*params->n_tables/128;
	}
	VectorU32 nsync_per_thread(params->n_threads);
#endif

	start_time = omp_get_wtime();
	omp_set_num_threads(params->n_threads); // split the windows across the threads
	#pragma omp parallel reduction(+:n_valid_windows, n_valid_hashes, n_bucket_entries, n_filtered)
	{
	    int tid = omp_get_thread_num();
	    int n_threads = omp_get_num_threads();
	    seq_t chunk_start = ((len - params->ref_window_size + 1) / n_threads)*tid;
	    seq_t chunk_end = ((len - params->ref_window_size + 1) / n_threads)*(tid + 1);
	    printf("Thread %d range: %u %u \n", tid, chunk_start, chunk_end);

#if EXTERNAL_MEM_INDEX
	    int sync_point = 1;
#endif
	    bool init_minhash = true;
	    for (seq_t pos = chunk_start; pos != chunk_end; pos++) { // for each window of the thread's chunk
	    	if((pos - chunk_start) % REPORT_WINDOW_PROC_GRANULARITY == 0 && (pos - chunk_start) != 0) {
				printf("Thread %d processed %u valid windows \n", tid, pos - chunk_start);
			}
#if EXTERNAL_MEM_INDEX
	    	// check if we should write to file
	    	if(file_nsync_points > 0) {
	    		int sync_chunk_size = (chunk_end - chunk_start + 1)/(file_nsync_points + 1);
	    		if(n_valid_windows == (uint32)sync_chunk_size*sync_point) {
	    			printf("Thread %d sync point: %u, n_valid_windows: %u \n", tid, pos, n_valid_windows);
	    			store_ref_idx_per_thread(tid, sync_point == 1, fastaFname, ref, params);
	    			sync_point++;
	    			nsync_per_thread[tid]++;
	    		}
	    	}
#endif
	    	// discard windows with low information content
	    	if(ignore_window_bitmask[pos]) {
	    		init_minhash = true;
	    		continue;
			}
	    	n_valid_windows++;

	    	// get the min-hash signature for the window
	    	VectorMinHash& minhashes = minhash_thread_vectors[tid]; // each thread indexes into its pre-allocated buffer
	    	minhash_matrix_t& rolling_minhash_matrix = minhash_matrices[tid];
	    	bool valid_hash;
	    	if(init_minhash == true) {
	    		valid_hash = minhash_rolling_init(seq.c_str(), pos, params->ref_window_size,
	    					rolling_minhash_matrix, ignore_kmer_bitmask, params,
							minhashes);
	    		init_minhash = false;
	    	} else {
	    		valid_hash = minhash_rolling(seq.c_str(), pos, params->ref_window_size,
	    					rolling_minhash_matrix, ignore_kmer_bitmask, params,
							minhashes);
	    	}

	    	if(!valid_hash) {
	    		continue;
	    	}
	    	n_valid_hashes++;

	    	for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
	    		minhash_t proj_hash = params->sketch_proj_hash_func.apply_vector(
	    				minhashes, params->sketch_proj_indices, t*params->sketch_proj_len);
	    		minhash_t bucket_hash = params->sketch_proj_hash_func.bucket_hash(proj_hash);
				buckets_t* buckets = &mutable_index.per_table_buckets[t];
	    		VectorSeqPos& bucket = buckets->per_thread_buckets_data_vectors[tid][bucket_hash];
	    		if(bucket.size() == 0) { // first item in this thread bucket
	    			bucket.resize(params->bucket_size);
	    		}
	    		// add to the bucket
	    		uint32 curr_size = buckets->per_thread_bucket_sizes[tid][bucket_hash];
	    		if(curr_size + 1 >= bucket.size()) {
	    			bucket.resize(curr_size + 100);
	    		}
				bool store_pos = true;
				if(curr_size > 0) {
					loc_t* epos = &bucket[curr_size-1];
					if(epos->len < MAX_LOC_LEN && (epos->pos + epos->len) == pos) {
						epos->len++;
						store_pos = false;
					} /*else if((epos->pos + epos->len) + params->bucket_entry_coverage >= pos) { // sampling
						epos->len += pos - (epos->pos + epos->len - 1);
						store_pos = false;
					}*/
				}
				if(store_pos) {
					loc_t new_loc;
					new_loc.pos = pos;
					new_loc.len = 1;
					new_loc.hash = proj_hash;
					bucket[curr_size] = new_loc;
					buckets->per_thread_bucket_sizes[tid][bucket_hash]++;
					n_bucket_entries++;
				} else {
					n_filtered++;
				}
	    	}
	    }
	}
	printf("Populated all the buckets. Time : %.2f sec\n", omp_get_wtime() - start_time);

	// 4. collect per thread results
	printf("Collecting thread buckets... \n");
	double start_coll_sort = omp_get_wtime();
	for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
		buckets_t* buckets = &mutable_index.per_table_buckets[t];
		buckets->n_entries = 0;
		for(uint32 b = 0; b < params->n_buckets; b++) {
			for(uint32 tid = 0; tid < params->n_threads; tid++) {
				VectorSeqPos& global_bucket = buckets->buckets_data_vectors[b];
				VectorSeqPos& thread_bucket = buckets->per_thread_buckets_data_vectors[tid][b];
				global_bucket.insert(global_bucket.end(), thread_bucket.begin(), thread_bucket.begin() + buckets->per_thread_bucket_sizes[tid][b]);
				VectorSeqPos().swap(buckets->per_thread_buckets_data_vectors[tid][b]);
				buckets->n_entries += buckets->per_thread_bucket_sizes[tid][b];
			}
		}
		std::vector<VectorU32>().swap(buckets->per_thread_bucket_sizes);
		std::vector<std::vector<VectorSeqPos>>().swap(buckets->per_thread_buckets_data_vectors);
	}
	printf("Collected all the buckets. Time : %.2f sec\n", omp_get_wtime() - start_coll_sort);

#if EXTERNAL_MEM_INDEX
	if(file_nsync_points > 0) { // read partial files from disk
		for(uint32 tid = 0; tid < params->n_threads; tid++) {
			load_ref_idx_per_thread(tid, nsync_per_thread[tid], fastaFname, ref, params);
		}
	}
#endif

	// 5. sort each bucket!
	printf("Sorting buckets... \n");
	double start_time_sort = omp_get_wtime();
	#pragma omp parallel for
	for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
		buckets_t* buckets = &mutable_index.per_table_buckets[t];
		for(uint32 b = 0; b < params->n_buckets; b++) {
			VectorSeqPos& bucket = buckets->buckets_data_vectors[b];
			std::sort(bucket.begin(), bucket.end(), comp_loc());
		}
	}
	printf("Total sort time : %.2f sec\n", omp_get_wtime() - start_time_sort);
	printf("Total number of valid reference windows: %u \n", n_valid_windows);
	printf("Total number of valid reference windows with valid hashes: %u \n", n_valid_hashes);
	printf("Total number of window bucket entries: %llu \n", n_bucket_entries);
	printf("Total number of window bucket entries filtered: %llu \n", n_filtered);
	printf("Total hashing time: %.2f sec\n", omp_get_wtime() - start_time);
}

void ref_t::load_index(const char* fastaFname, const index_params_t* params) {
	double start_time = omp_get_wtime();
	fasta2ref(fastaFname, *this);
	load_freq_kmers(fastaFname, high_freq_kmer_bitmap, high_freq_kmer_trie, params->max_count);
	if(params->load_mhi) {
		load_ref_idx(fastaFname, *this, params);
	}
	printf("Index loading time: %.2f sec\n", (float)(omp_get_wtime() - start_time));
}

void ref_t::store_index(const char* fastaFname, const index_params_t* params) {
	double start_time = omp_get_wtime();
	store_ref_idx(fastaFname, *this, params);
	printf("Index storing time: %.2f sec\n", (float)(omp_get_wtime() - start_time));
}

void ref_t::mark_freq_kmers(const index_params_t* params) {
	clock_t t = clock();
	ignore_kmer_bitmask.resize(len - params->k);
	#pragma omp parallel for
	for(seq_t i = 0; i < len - params->k + 1; i++) { // for each window of the genome
		uint32_t packed_kmer;
		if(pack_32(&seq.c_str()[i], params->k, &packed_kmer) < 0) {
			ignore_kmer_bitmask[i] = 1; // contains ambiguous bases
			continue;
		}
		if(high_freq_kmer_bitmap[packed_kmer]) {
			ignore_kmer_bitmask[i] = 1;
		}
	}
	printf("Done marking frequent kmers time: %.2f sec \n", (float) (clock() - t)/CLOCKS_PER_SEC);
}

// checks if the given sequence is informative or not
// e.g. non-informative seq: same character is repeated throughout the seq (NN...N)
#define AMBIG_BASE_FRAC 50
#define LOW_BASE_FRAC 100
int is_inform_ref_window(const char* seq, const uint32_t len, const index_params_t* params) {
	uint32 base_counts[5] = { 0 };
	for(uint32 i = 0; i < len; i++) {
		base_counts[(int) seq[i]]++;
	}
	if(base_counts[4] > params->ref_window_size/AMBIG_BASE_FRAC) { // N ambiguous bases
		return 0;
	}
	uint32 n_empty = 0;
	uint32 n_low = 0;
	for(uint32 i = 0; i < 4; i++) {
		if(base_counts[i] == 0) {
			n_empty++;
		}
		if(base_counts[i] < params->ref_window_size/LOW_BASE_FRAC) {
			n_low++;
		}
	}
	if(n_empty > 1 || n_low > 1) { // repetitions of 2 or 1 base
		return 0;
	}
	return 1;
}

void ref_t::mark_windows_to_discard(const index_params_t* params) {
	ignore_window_bitmask.resize(len - params->ref_window_size + 1);
	#pragma omp parallel for
	for(seq_t pos = 0; pos < len - params->ref_window_size + 1; pos++) { // for each window of the genome
		if(!is_inform_ref_window(&seq.c_str()[pos], params->ref_window_size, params)) {
			ignore_window_bitmask[pos] = 1; // discard windows with low information content
		}
	}
}
