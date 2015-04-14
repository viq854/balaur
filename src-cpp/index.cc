#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <time.h>

#include <omp.h>
#include "types.h"
#include "index.h"
#include "lsh.h"
#include "io.h"
#include "hash.h"


// checks if the given sequence is informative or not
// e.g. non-informative seq: same character is repeated throughout the seq (NN...N)
int is_inform_ref_window(const char* seq, const uint32_t len, const index_params_t* params) {
	uint32 base_counts[5] = { 0 };
	for(uint32 i = 0; i < len; i++) {
		base_counts[(int) seq[i]]++;
	}
	if(base_counts[4] > params->ref_window_size/50) { // N ambiguous bases
		return 0;
	}
	uint32 n_empty = 0;
	for(uint32 i = 0; i < 4; i++) {
		if(base_counts[i] == 0) {
			n_empty++;
		}
	}
	if(n_empty > 1) { // repetitions of 2 or 1 base
		return 0;
	}

	uint32 n_low = 0;
	for(uint32 i = 0; i < 4; i++) {
		if(base_counts[i] < params->ref_window_size/100) {
			n_low++;
		}
	}

	if(n_low > 1) { // almost repetitions of 2 or 1 base
		return 0;
	}

	return 1;
}

// compute and store the frequency of each kmer in the given sequence
void compute_kmer_counts(const char* seq, const seq_t seq_len, const index_params_t* params,
		MapKmerCounts& hist) {
	for(seq_t j = 0; j <= (seq_len - params->k); j++) {
		uint32_t kmer;
		if(pack_32(&seq[j], params->k, &kmer) < 0) {
			continue;
		}
		hist[kmer]++;
	}
}

void mark_freq_kmers(ref_t& ref, const index_params_t* params) {
	clock_t t = clock();
	ref.ignore_kmer_bitmask.resize(ref.len - params->k);
	omp_set_num_threads(params->n_threads); // split the windows across the threads
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int n_threads = omp_get_num_threads();
		seq_t chunk_start = tid*(ref.len - params->k + 1) / n_threads;
		seq_t chunk_end = (tid + 1)*(ref.len - params->k + 1) / n_threads;

		marisa::Agent agent;
		for(seq_t i = chunk_start; i < chunk_end; i++) {
			for (uint32 k = 0; k < params->k; k++) {
				if(ref.seq.c_str()[i+k] == BASE_IGNORE) {
					ref.ignore_kmer_bitmask[i] = 1; // contains ambiguous bases
					break;
				}
			}
			agent.set_query(&ref.seq.c_str()[i], params->k);
			if(ref.high_freq_kmer_trie.lookup(agent)) {
				ref.ignore_kmer_bitmask[i] = 1;
			}
		}
	}
	printf("Done marking frequent kmers time: %.2f sec \n", (float) (clock() - t)/CLOCKS_PER_SEC);
}

void mark_windows_to_discard(ref_t& ref, const index_params_t* params) {
	ref.ignore_window_bitmask.resize(ref.len - params->ref_window_size + 1);
	#pragma omp parallel for
	for(seq_t pos = 0; pos < ref.len - params->ref_window_size + 1; pos++) { // for each window of the genome
		if(!is_inform_ref_window(&ref.seq.c_str()[pos], params->ref_window_size, params)) {
			ref.ignore_window_bitmask[pos] = 1; // discard windows with low information content
		}
	}
}

// --- Minhash ---

// index reference:
// - load fasta
// - generate valid windows (sliding window)
// - compute the hash of each window
// - sort

#define MAX_NTABLES_NO_DISK 1024

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
	//compute_kmer_counts(ref.seq.c_str(), ref.seq.size(), params, ref.kmer_hist);
	//store_kmer_hist(fastaFname, ref.kmer_hist);
	//ref.kmer_hist = MapKmerCounts(); // free memory
	load_freq_kmers(fastaFname, ref.high_freq_kmer_trie, params->max_count);
	mark_freq_kmers(ref, params);

	printf("Loading valid windows mask... \n");
	if(!load_valid_window_mask(fastaFname, ref, params)) {
		mark_windows_to_discard(ref, params);
		store_valid_window_mask(fastaFname, ref, params);
	}
	printf("Total window/kmer pre-processing time: %.2f sec\n", omp_get_wtime() - start_time);

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
	std::vector<VectorMinHash> minhash_thread_vectors(params->n_threads);
	for(uint32 i = 0; i < params->n_threads; i++) {
		minhash_thread_vectors[i].resize(params->h);
		std::fill(minhash_thread_vectors[i].begin(), minhash_thread_vectors[i].end(), UINT_MAX);
	}
	std::vector<CyclicHash*> hasher_thread_vectors(params->n_threads);
	for(uint32 i = 0; i < params->n_threads; i++) {
		hasher_thread_vectors[i] = new CyclicHash(params->k, 32);
	}
	std::vector<minhash_matrix_t> minhash_matrices(params->n_threads);

	// 3. hash each valid window
	printf("Hashing reference windows... \n");
	uint32 n_valid_windows = 0;
	uint32 n_valid_hashes = 0;
	uint64 n_bucket_entries = 0;
	uint64 n_filtered = 0;
	uint64 n_dropped = 0;

	int file_nsync_points = 0;
	if(params->n_tables > MAX_NTABLES_NO_DISK) {
		file_nsync_points = 42*params->n_tables/128; //params->n_tables / MAX_NTABLES_NO_DISK - 1;
	}
	VectorU32 nsync_per_thread(params->n_threads);

	start_time = omp_get_wtime();
	omp_set_num_threads(params->n_threads); // split the windows across the threads
	#pragma omp parallel reduction(+:n_valid_windows, n_valid_hashes, n_bucket_entries, n_filtered, n_dropped)
	{
	    int tid = omp_get_thread_num();
	    int n_threads = omp_get_num_threads();
	    seq_t chunk_start = ((ref.len - params->ref_window_size + 1) / n_threads)*tid;
	    seq_t chunk_end = ((ref.len - params->ref_window_size + 1) / n_threads)*(tid + 1);
	    printf("Thread %d range: %u %u \n", tid, chunk_start, chunk_end);

	    int sync_point = 1;
	    bool init_minhash = true;
	    for (seq_t pos = chunk_start; pos != chunk_end; pos++) { // for each window of the thread's chunk

	    	// check if we should write to file
	    	if(file_nsync_points > 0) {
	    		int sync_chunk_size = (chunk_end - chunk_start + 1)/(file_nsync_points + 1);
	    		if(n_valid_windows == sync_chunk_size*sync_point) {
	    			printf("Thread %d sync point: %u, n_valid_windows: %u \n", tid, pos, n_valid_windows);
	    			store_ref_idx_per_thread(tid, sync_point == 1, fastaFname, ref, params);
	    			sync_point++;
	    			nsync_per_thread[tid]++;
	    		}
	    	}

	    	if(ref.ignore_window_bitmask[pos]) { // discard windows with low information content
	    		init_minhash = true;
	    		continue;
			}
	    	n_valid_windows++;
	    	if((pos - chunk_start) % 2000000 == 0 && (pos - chunk_start) != 0) {
	    		printf("Thread %d processed %u valid windows \n", tid, pos - chunk_start);
	    	}

	    	// get the min-hash signature for the window
	    	VectorMinHash& minhashes = minhash_thread_vectors[tid]; // each thread indexes into its pre-allocated buffer
	    	minhash_matrix_t& rolling_minhash_matrix = minhash_matrices[tid];
	    	bool valid_hash;
	    	if(init_minhash == true) { //if(pos == chunk_start) {
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
	    	if(!valid_hash) {
	    		continue;
	    	}
	    	n_valid_hashes++;

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
	    			bucket.resize(curr_size + 200);
	    		}
				bool store_pos = true;
				if(curr_size > 0) {
					loc_t* epos = &bucket[curr_size-1];
					if((epos->pos + epos->len) == pos) {
						epos->len++;
						store_pos = false;
					} /*else if((epos->pos + epos->len) + params->bucket_entry_coverage >= pos) {
						epos->len += pos - (epos->pos + epos->len - 1);
						store_pos = false;
					}*/
				}
				if(store_pos) {
					loc_t new_loc;
					new_loc.pos = pos;
					new_loc.len = 1;
					bucket[curr_size] = new_loc;
					buckets->per_thread_bucket_sizes[tid][bucket_index]++;
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
	if(file_nsync_points > 0) { // read partial files from disk
		for(uint32 tid = 0; tid < params->n_threads; tid++) {
			load_ref_idx_per_thread(tid, nsync_per_thread[tid], fastaFname, ref, params);
		}
	}

	printf("Collected all the buckets. Time : %.2f sec\n", omp_get_wtime() - start_coll_sort);

	// 5. sort each bucket!
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
	printf("Total number of window bucket entries dropped: %llu \n", n_dropped);
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

	// 3. load the frequency of each kmer and collect high-frequency kmers
	double start_time = omp_get_wtime();
	load_freq_kmers(fastaFname, ref.high_freq_kmer_trie, params->max_count);
	//mark_freq_kmers(ref, params);
	double end_time = omp_get_wtime();
	printf("Total kmer pre-processing time: %.2f sec\n", end_time - start_time);

	if(!load_kmer2_hashes(fastaFname, ref, params)) {
		compute_store_kmer2_hashes(fastaFname, ref, params);
	}
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
			//compute_kmer_counts(r.seq.c_str(), r.len, params, reads.kmer_hist);
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

		r->minhashes_rc.resize(params->h);
		r->valid_minhash_rc = minhash(r->rc.c_str(), 0, r->len,
				ref.high_freq_kmer_trie,
				ref.ignore_kmer_bitmask,
				marisa::Trie(), params,
				params->kmer_hasher, false,
				r->minhashes_rc);
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
