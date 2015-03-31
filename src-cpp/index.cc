#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <time.h>
#include "types.h"
#include "index.h"
#include "lsh.h"
#include "io.h"
#include "hash.h"

// compute and store the frequency of each kmer in the given sequence
void compute_kmer_histogram(const char* seq, const seq_t seq_len, const index_params_t* params, MapKmerCounts& hist) {
	for(int j = 0; j <= (seq_len - params->k); j++) {
		uint32_t kmer;
		if(pack_32(&seq[j], params->k, &kmer) < 0) {
			continue;
		}
		hist[kmer]++;
	}
}

// returns true if the kmer is uninformative
bool filter_kmer(const std::string& seq, const seq_t i, const marisa::Trie& high_freq_kmer_trie, const index_params_t* params) {
	for (uint32 k = 0; k < params->k; k++) {
		if(seq.c_str()[i+k] == BASE_IGNORE) {
			return true; // contains ambiguous bases
		}
	}
	marisa::Agent query_agent;
	query_agent.set_query(&seq.c_str()[i], params->k);
	if(high_freq_kmer_trie.lookup(query_agent)) {
		return true; // high-frequency kmer
	}
	return false;
}

// populates the forward and reverse complement filtered kmer bitmasks
void mark_kmers_to_discard(ref_t& ref, const index_params_t* params) {
	ref.ignore_kmer_bitmask.resize(ref.len - params->k);
	if(INDEX_READS_REF) ref.ignore_kmer_bitmask_RC.resize(ref.len - params->k);
	#pragma omp parallel for
	for(seq_t i = 0; i < ref.len - params->k + 1; i++) {
		if(filter_kmer(ref.seq, i, ref.high_freq_kmer_trie, params)) {
			ref.ignore_kmer_bitmask[i] = true;
			if(INDEX_READS_REF) ref.ignore_kmer_bitmask_RC[ref.len - i - params->k] = true;
		}
	}
}

// checks if the given sequence is informative or not
// e.g. non-informative seq: same character is repeated throughout the seq (NN...N)
bool filter_window(const char* seq, const uint32_t len) {
	uint32 base_counts[5] = { 0 };
	for(uint32 i = 0; i < len; i++) {
		base_counts[(int) seq[i]]++;
	}
	if(base_counts[4] > 10) { // N ambiguous bases
		return true;
	}
	uint32 n_empty = 0;
	for(uint32 i = 0; i < 4; i++) {
		if(base_counts[i] == 0) {
			n_empty++;
		}
	}
	if(n_empty > 1) { // repetitions of 2 or 1 base
		return true;
	}

	return false;
}

void mark_windows_to_discard(ref_t& ref, const index_params_t* params) {
	ref.ignore_window_bitmask.resize(ref.len - params->ref_window_size + 1);
	if(INDEX_READS_REF) ref.ignore_window_bitmask_RC.resize(ref.len - params->ref_window_size + 1);
	#pragma omp parallel for
	for(seq_t pos = 0; pos < ref.len - params->ref_window_size + 1; pos++) { // for each window of the genome
		if(filter_window(&ref.seq.c_str()[pos], params->ref_window_size)) {
			ref.ignore_window_bitmask[pos] = true; // discard windows with low information content
			if(INDEX_READS_REF) ref.ignore_window_bitmask_RC[ref.len - pos - params->ref_window_size] = true;
		}
	}
}

// --- Minhash ---

// index reference:
// - load fasta
// - generate valid windows (sliding window)
// - compute the hash of each window
// - sort

int add_window_to_buckets(const seq_t window_pos, const int tid, const VectorMinHash& minhashes, const index_params_t* params, ref_t& ref, const uint16_t rc) {
	int n_new_entries = 0;
	for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
		minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(
				minhashes,
				params->sketch_proj_indices,
				t*params->sketch_proj_len);

		buckets_t* buckets = &ref.hash_tables[t];
		uint32 bucket_index = buckets->per_thread_bucket_indices[tid][bucket_hash];
		if(bucket_index == buckets->n_buckets) { // this is the first entry in the bucket
			bucket_index = buckets->per_thread_next_free_bucket_index[tid];
			buckets->per_thread_next_free_bucket_index[tid]++;
			buckets->per_thread_bucket_indices[tid][bucket_hash] = bucket_index;
		}

		VectorSeqPos& bucket = buckets->per_thread_buckets_data_vectors[tid][bucket_index];
		if(bucket.size() == 0) { // first item in this thread bucket
			bucket.resize(params->bucket_size);
		}
		// add to the bucket
		uint32 curr_size = buckets->per_thread_bucket_sizes[tid][bucket_index];
		if(curr_size + 1 >= bucket.size()) {
			bucket.resize(curr_size + 500);
		}
		if(curr_size > 0) {
			loc_t* epos = &bucket[curr_size-1];
			if((epos->pos + epos->len) == window_pos) {
				epos->len++; // extend the existing contig
				continue;
			}
		}
		// add a new entry
		loc_t new_loc;
		new_loc.pos = window_pos;
		new_loc.rc = rc;
		new_loc.len = 1;
		bucket[curr_size] = new_loc;
		buckets->per_thread_bucket_sizes[tid][bucket_index]++;
		n_new_entries++;
	}
	return n_new_entries;
}

void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref) {
	printf("**** SRX Reference Indexing ****\n");
		
	// 1. load the reference
	printf("Loading FASTA file %s... \n", fastaFname);
	clock_t t = clock();
	fasta2ref(fastaFname, ref);
	printf("Reference loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 2. load the frequency of each kmer and collect high-frequency kmers
	printf("Loading frequent kmers... \n");
	double start_time = omp_get_wtime();
	//compute_kmer_histogram(ref.seq.c_str(), ref.seq.size(), params, ref.kmer_hist);
	//store_kmer_hist(fastaFname, ref.kmer_hist);
	//ref.kmer_hist = MapKmerCounts(); // free memory
	load_freq_kmers(fastaFname, ref.high_freq_kmer_trie, ref.high_freq_kmer_trie_RC, params->max_count);
	mark_kmers_to_discard(ref, params);
	printf("Total kmer pre-processing time: %.2f sec\n", omp_get_wtime() - start_time);

	//mark_windows_to_discard(ref, params);
	//store_valid_window_mask(fastaFname, ref);
	load_valid_window_mask(fastaFname, ref, params);

	// initialize the hash tables
	ref.hash_tables.resize(params->n_tables);
	for(uint32 t = 0; t < params->n_tables; t++) {
		buckets_t* buckets = &ref.hash_tables[t];
		buckets->n_buckets = pow(2, params->n_buckets_pow2);

		// initialize global buckets
		buckets->next_free_bucket_index = 0;
		buckets->bucket_indices.resize(buckets->n_buckets, buckets->n_buckets);
		buckets->buckets_data_vectors.resize(buckets->n_buckets);

		// initialize per thread buckets
		buckets->per_thread_next_free_bucket_index.resize(params->n_threads, 0);
		buckets->per_thread_bucket_indices.resize(params->n_threads);
		buckets->per_thread_buckets_data_vectors.resize(params->n_threads);
		buckets->per_thread_bucket_sizes.resize(params->n_threads);
		for(uint32 tid = 0; tid < params->n_threads; tid++) {
			buckets->per_thread_bucket_indices[tid].resize(buckets->n_buckets, buckets->n_buckets);
			buckets->per_thread_buckets_data_vectors[tid].resize(buckets->n_buckets);
			buckets->per_thread_bucket_sizes[tid].resize(buckets->n_buckets, 0);
		}
	}

	// initialize additional per-thread storage
	std::vector<minhash_matrix_t> minhash_matrices(params->n_threads);
	std::vector<VectorMinHash> minhash_thread_vectors(params->n_threads);
	for(uint32 i = 0; i < params->n_threads; i++) {
		minhash_thread_vectors[i].resize(params->h);
		std::fill(minhash_thread_vectors[i].begin(), minhash_thread_vectors[i].end(), UINT_MAX);
	}
	std::vector<CyclicHash*> hasher_thread_vectors(params->n_threads);
	for(uint32 i = 0; i < params->n_threads; i++) {
		hasher_thread_vectors[i] = new CyclicHash(params->k, 32);
	}

	// 3. hash each valid window
	printf("Hashing reference windows... \n");
	uint32 n_valid_windows = 0;
	uint32 n_valid_hashes = 0;
	uint64 n_bucket_entries = 0;
	uint64 n_filtered = 0;

	int file_nsync_points = 0;
	if(DISK_SYNC_PARTIAL_TABLES) {
		file_nsync_points = 52*params->n_tables/128;
	}
	VectorU32 nsync_per_thread(params->n_threads);

	start_time = omp_get_wtime();
	omp_set_num_threads(params->n_threads); // split the windows across the threads
	#pragma omp parallel reduction(+:n_valid_windows, n_valid_hashes, n_bucket_entries, n_filtered)
	{
	    int tid = omp_get_thread_num();
	    int n_threads = omp_get_num_threads();
	    seq_t chunk_start = ((ref.len - params->ref_window_size + 1) / n_threads)*tid;
	    seq_t chunk_end = ((ref.len - params->ref_window_size + 1) / n_threads)*(tid + 1);
	    printf("Thread %d windows: %u - %u \n", tid, chunk_start, chunk_end);

	    VectorMinHash& minhashes = minhash_thread_vectors[tid]; // each thread indexes into its pre-allocated buffer
	    minhash_matrix_t& rolling_minhash_matrix = minhash_matrices[tid];
	    int sync_point = 1;
	    bool init_minhash = true;
	    for(seq_t pos = chunk_start; pos != chunk_end; pos++) { // for each window of the thread's chunk
	    	if((pos - chunk_start) % 2000000 == 0 && (pos - chunk_start) != 0) {
				printf("[PROGRESS] Thread %d processed %u valid windows \n", tid, pos - chunk_start);
			}

	    	if(DISK_SYNC_PARTIAL_TABLES && file_nsync_points > 0) { // check if we should write partial results to disk
	    		int sync_chunk_size = (chunk_end - chunk_start + 1)/(file_nsync_points + 1);
	    		if(n_valid_windows == sync_chunk_size*sync_point) {
	    			printf("[DISK] Thread %d sync point: %u, n_valid_windows: %u \n", tid, pos, n_valid_windows);
	    			store_ref_idx_per_thread(tid, sync_point == 1, fastaFname, ref, params);
	    			sync_point++;
	    			nsync_per_thread[tid]++;
	    		}
	    	}

	    	// FORWARD
	    	if(ref.ignore_window_bitmask[pos]) { // discard windows with low information content
	    		init_minhash = true;
			} else {
				n_valid_windows++;
				bool valid_hash;
				if(init_minhash == true) {
					valid_hash = minhash_rolling_init(ref.seq.c_str(), pos, params->ref_window_size,
								rolling_minhash_matrix,
								ref.ignore_kmer_bitmask, params,
								hasher_thread_vectors[0],
								minhashes);
					init_minhash = false;
				} else {
					valid_hash = minhash_rolling(ref.seq.c_str(), pos, params->ref_window_size,
								rolling_minhash_matrix,
								ref.ignore_kmer_bitmask, params,
								hasher_thread_vectors[0],
								minhashes);
				}
				if(valid_hash) {
					n_valid_hashes++;
					int n_new_entries = add_window_to_buckets(pos, tid, minhashes, params, ref, 0);
					n_bucket_entries += n_new_entries;
					n_filtered += params->n_tables - n_new_entries;
				}
			}
	    }
	    // REVERSE COMPLEMENT
	    if(INDEX_READS_REF) {
	    	std::fill(minhash_thread_vectors[tid].begin(), minhash_thread_vectors[tid].end(), UINT_MAX);
	    	init_minhash = true;
	    	for(seq_t pos = chunk_start; pos != chunk_end; pos++) { // for each window of the thread's RC chunk
	    		if((pos - chunk_start) % 2000000 == 0 && (pos - chunk_start) != 0) {
	    			printf("[PROGRESS] Thread %d processed %u valid RC windows \n", tid, pos - chunk_start);
	    		}

	    		if(DISK_SYNC_PARTIAL_TABLES && file_nsync_points > 0) { // check if we should write partial results to disk
	    			int sync_chunk_size = (chunk_end - chunk_start + 1)/(file_nsync_points + 1);
	    			if(n_valid_windows == sync_chunk_size*sync_point) {
	    				printf("[DISK] Thread %d sync point: %u, n_valid_windows: %u \n", tid, pos, n_valid_windows);
	    				store_ref_idx_per_thread(tid, sync_point == 1, fastaFname, ref, params);
	    				sync_point++;
	    				nsync_per_thread[tid]++;
	    			}
	    		}
	    		if(ref.ignore_window_bitmask_RC[pos]) { // discard windows with low information content
	    			init_minhash = true;
	    		} else {
	    			n_valid_windows++;
					bool valid_hash;
					if(init_minhash == true) {
						valid_hash = minhash_rolling_init(ref.seq_RC.c_str(), pos, params->ref_window_size,
									rolling_minhash_matrix,
									ref.ignore_kmer_bitmask_RC, params,
									hasher_thread_vectors[0],
									minhashes);
						init_minhash = false;
					} else {
						valid_hash = minhash_rolling(ref.seq_RC.c_str(), pos, params->ref_window_size,
									rolling_minhash_matrix,
									ref.ignore_kmer_bitmask_RC, params,
									hasher_thread_vectors[0],
									minhashes);
					}
					if(valid_hash) {
						add_window_to_buckets(pos, tid, minhashes, params, ref, 1);
					}
	    		}
	    	}
	    }
	}
	printf("Populated all the buckets. Time : %.2f sec\n", omp_get_wtime() - start_time);

	// 4. collect per thread results
	printf("Collecting thread buckets... \n");
	double start_coll_sort = omp_get_wtime();
	for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
		buckets_t* buckets = &ref.hash_tables[t];
		for(uint32 b = 0; b < buckets->n_buckets; b++) {
			for(uint32 tid = 0; tid < params->n_threads; tid++) {
				uint32 thread_bucket_index = buckets->per_thread_bucket_indices[tid][b];
				if(thread_bucket_index == buckets->n_buckets) continue;
				uint32 global_bucket_index = buckets->bucket_indices[b];
				if(global_bucket_index == buckets->n_buckets) {
					global_bucket_index = buckets->next_free_bucket_index;
					buckets->next_free_bucket_index++;
					buckets->bucket_indices[b] = global_bucket_index;
				}
				VectorSeqPos& global_bucket = buckets->buckets_data_vectors[global_bucket_index];
				VectorSeqPos& thread_bucket = buckets->per_thread_buckets_data_vectors[tid][thread_bucket_index];
				global_bucket.insert(global_bucket.end(), thread_bucket.begin(), thread_bucket.begin() + buckets->per_thread_bucket_sizes[tid][thread_bucket_index]);
				VectorSeqPos().swap(buckets->per_thread_buckets_data_vectors[tid][thread_bucket_index]);
			}
		}
		//VectorPerThreadIndices().swap(buckets->per_thread_bucket_indices);
		VectorPerThreadSizes().swap(buckets->per_thread_bucket_sizes);
		VectorU32().swap(buckets->per_thread_next_free_bucket_index);
		VectorPerThreadBuckets().swap(buckets->per_thread_buckets_data_vectors);
	}
	if(DISK_SYNC_PARTIAL_TABLES && file_nsync_points > 0) { // read partial files from disk
		for(uint32 tid = 0; tid < params->n_threads; tid++) {
			load_ref_idx_per_thread(tid, nsync_per_thread[tid], fastaFname, ref, params);
		}
	}
	printf("Collected all the buckets. Time : %.2f sec\n", omp_get_wtime() - start_coll_sort);

	// 5. sort each bucket
	printf("Sorting buckets... \n");
	double start_time_sort = omp_get_wtime();
	#pragma omp parallel for
	for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
		buckets_t* buckets = &ref.hash_tables[t];
		for(uint32 b = 0; b < buckets->next_free_bucket_index; b++) {
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

void load_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref) {
	printf("**** SRX Reference Index Loading ****\n");

	// 1. load the reference
	printf("Loading FASTA file %s... \n", fastaFname);
	clock_t t = clock();
	fasta2ref(fastaFname, ref);
	printf("Reference loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 2. load the index
	printf("Loading reference index for reference file %s... \n", fastaFname);
	t = clock();
	load_ref_idx(fastaFname, ref, params);
	printf("Reference index loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 3. load the high-frequency kmers
	printf("Loading high-frequency kmers trie... \n");
	t = clock();
	load_freq_kmers(fastaFname, ref.high_freq_kmer_trie, ref.high_freq_kmer_trie_RC, params->max_count);
	printf("Total kmer trie load time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}


void store_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref) {
	printf("**** SRX Reference Index Loading ****\n");

	// store the reference buckets
	printf("Storing the reference index for reference file %s... \n", fastaFname);
	clock_t t = clock();
	store_ref_idx(fastaFname, ref, params);
	printf("Reference storing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}

void index_reads_lsh(const char* readsFname, ref_t& ref, index_params_t* params, reads_t& reads) {
	printf("**** SRX Read Indexing ****\n");

	// 1. load the reads (TODO: batch mode)
	printf("Loading FATQ reads...\n");
	clock_t t = clock();
	fastq2reads(readsFname, reads);
	printf("Total read loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 2. compute the frequency of each kmer in the read set
	t = clock();
	if(params->min_count != 0) {
		/*for(uint32 i = 0; i < reads.reads.size(); i++) { // add the contribution of each read
			read_t r = reads.reads[i];
			//compute_kmer_histogram(r.seq.c_str(), r.len, params, reads.kmer_hist);
		}
		find_low_freq_kmers(reads.kmer_hist, reads.low_freq_kmer_hist, params);
		reads.kmer_hist = MapKmerCounts(); // free memory
		*/
	}
	printf("Total kmer pre-processing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 3. compute the min-hash signature of each read
	printf("Hashing reads...\n");
	double start_time = omp_get_wtime();
	omp_set_num_threads(params->n_threads);
	#pragma omp parallel for
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		r->minhashes.resize(params->h);
		r->valid_minhash = minhash(r->seq.c_str(), 0, r->len,
				ref.high_freq_kmer_trie,
				ref.ignore_kmer_bitmask,
				marisa::Trie(), params,
				params->kmer_hasher, false,
				r->minhashes);

		if(INDEX_READS_RC) {
			r->minhashes_rc.resize(params->h);
			r->valid_minhash_rc = minhash(r->rc.c_str(), 0, r->len,
					ref.high_freq_kmer_trie,
					ref.ignore_kmer_bitmask,
					marisa::Trie(), params,
					params->kmer_hasher, false,
					r->minhashes_rc);
		} else {
			r->valid_minhash_rc = 0;
		}
	}
	double end_time = omp_get_wtime();
	printf("Total read hashing time: %.2f sec\n", end_time - start_time);
	// TODO: free memory
}

// --- Sampling ---

// create hash table i 
// - hash each window
// - sort
//void index_ref_table_i(ref_t& ref, const index_params_t* params, const seq_t i) {
//	// hash each window using sampling ids i
//	clock_t t = clock();
//	for(seq_t j = 0; j < ref.windows_by_pos.size(); j++) {
//		ref.windows_by_pos[j].simhash = sampling_hash(ref.seq.c_str(), ref.windows_by_pos[j].pos, i, params);
//	}
//	printf("Total hash table %llu computation time: %.2f sec\n", (uint64) i, (float)(clock() - t) / CLOCKS_PER_SEC);
//
//	// sort the hashes
//	t = clock();
//	sort_windows_hash(ref);
//	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
//	printf("Total number of distinct window hashes = %llu \n", (uint64) get_num_distinct(ref));
//}
//
//void index_reads_table_i(reads_t& reads, const index_params_t* params, const seq_t i) {
//	// hash each read using sampling ids i
//	clock_t t = clock();
//	for(uint32 j = 0; j < reads.reads.size(); j++) {
//		reads.reads[j].simhash = sampling_hash(reads.reads[j].seq.c_str(), 0, i, params);
//	}
//	printf("Total read hash table %llu computation time: %.2f sec\n", (uint64) i, (float)(clock() - t) / CLOCKS_PER_SEC);
//}
