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
int is_inform_ref_window(const char* seq, const uint32_t len) {
	uint32 base_counts[5] = { 0 };
	for(uint32 i = 0; i < len; i++) {
		base_counts[(int) seq[i]]++;
	}
	if(base_counts[4] > 10) { // N ambiguous bases
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
	ref.ignore_window_bitmask.resize(ref.len- params->ref_window_size + 1);
	#pragma omp parallel for
	for(seq_t pos = 0; pos < ref.len - params->ref_window_size + 1; pos++) { // for each window of the genome
		if(!is_inform_ref_window(&ref.seq.c_str()[pos], params->ref_window_size)) {
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
void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref) {
	printf("**** SRX Reference Indexing ****\n");
		
	// 1. load the reference
	printf("Loading FASTA file %s... \n", fastaFname);
	clock_t t = clock();
	fasta2ref(fastaFname, ref);
	printf("Reference loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 2. load the frequency of each kmer and collect high-frequency kmers
	double start_time = omp_get_wtime();
	//compute_kmer_counts(ref.seq.c_str(), ref.seq.size(), params, ref.kmer_hist);
	//store_kmer_hist(fastaFname, ref.kmer_hist);
	//ref.kmer_hist = MapKmerCounts(); // free memory
	load_freq_kmers(fastaFname, ref.high_freq_kmer_trie, params->max_count);
	mark_freq_kmers(ref, params);
	double end_time = omp_get_wtime();
	printf("Total kmer pre-processing time: %.2f sec\n", end_time - start_time);

	// 3. hash each valid window

	// initialize the hash tables
	ref.hash_tables.resize(params->n_tables);
	for(uint32 t = 0; t < params->n_tables; t++) {
		buckets_t* buckets = &ref.hash_tables[t];
		buckets->n_buckets = pow(2, params->n_buckets_pow2);
		buckets->next_free_bucket_index = 0;
		buckets->bucket_indices.resize(buckets->n_buckets, buckets->n_buckets);
		buckets->bucket_index_locks.resize(buckets->n_buckets);
		for(uint32 l = 0; l < buckets->n_buckets; l++) {
			omp_init_lock(&buckets->bucket_index_locks[l]);
		}
		omp_init_lock(&buckets->lock);
		buckets->buckets_data_vectors.resize(buckets->n_buckets);
		buckets->bucket_sizes.resize(buckets->n_buckets, 0);

		// per thread buckets
		buckets->per_thread_buckets_data_vectors.resize(buckets->n_buckets);
		buckets->per_thread_bucket_sizes.resize(buckets->n_buckets);
		for(uint32 i = 0; i < buckets->n_buckets; i++) {
			buckets->per_thread_buckets_data_vectors[i].resize(params->n_threads);
		}
		for(uint32 i = 0; i < buckets->n_buckets; i++) {
			buckets->per_thread_bucket_sizes[i].resize(params->n_threads, 0);
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

	printf("Hashing reference windows... \n");
	uint32 n_valid_windows = 0;
	uint32 n_valid_hashes = 0;
	uint64 n_bucket_entries = 0;
	uint64 n_filtered = 0;
	uint64 n_dropped = 0;
	start_time = omp_get_wtime();
	omp_set_num_threads(params->n_threads); // split the windows across the threads
	#pragma omp parallel reduction(+:n_valid_windows, n_valid_hashes, n_bucket_entries, n_filtered, n_dropped)
	{
	    int tid = omp_get_thread_num();
	    int n_threads = omp_get_num_threads();
	    seq_t chunk_start = ((ref.len - params->ref_window_size + 1) / n_threads)*tid;
	    seq_t chunk_end = ((ref.len - params->ref_window_size + 1) / n_threads)*(tid + 1);
	    printf("Thread %d range: %u %u \n", tid, chunk_start, chunk_end);

	    for (seq_t pos = chunk_start; pos != chunk_end; pos++) { // for each window of the thread's chunk
	    	/*if(ref.ignore_window_bitmask[pos]) { // discard windows with low information content
				continue;
			}*/
	    	n_valid_windows++;
	    	if((pos - chunk_start) % 2000000 == 0 && (pos - chunk_start) != 0) {
	    		printf("Thread %d processed %u valid windows \n", tid, pos - chunk_start);
	    	}

	    	// get the min-hash signature for the window
	    	VectorMinHash& minhashes = minhash_thread_vectors[tid]; // each thread indexes into its pre-allocated buffer
	    	minhash_matrix_t& rolling_minhash_matrix = minhash_matrices[tid];
	    	bool valid_hash;
	    	if(pos == chunk_start) {
	    		valid_hash = minhash_rolling_init(ref.seq.c_str(), pos, params->ref_window_size,
	    					rolling_minhash_matrix,
							ref.ignore_kmer_bitmask, params,
							hasher_thread_vectors[0],
							minhashes);
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
	    		omp_set_lock(&buckets->bucket_index_locks[bucket_hash]);
	    		uint32 bucket_index = buckets->bucket_indices[bucket_hash];
	    		if(bucket_index == buckets->n_buckets) { // this is the first entry in the global bucket
	    			omp_set_lock(&buckets->lock);
	    			bucket_index = buckets->next_free_bucket_index;
	    			buckets->next_free_bucket_index++;
	    			buckets->bucket_indices[bucket_hash] = bucket_index;
	    			omp_unset_lock(&buckets->lock);
	    		}
	    		omp_unset_lock(&buckets->bucket_index_locks[bucket_hash]);
	    		VectorSeqPos& bucket = buckets->per_thread_buckets_data_vectors[bucket_index][tid];
	    		if(bucket.size() == 0) {
	    			// first item in this thread bucket
	    			bucket.resize(params->bucket_size);
	    		}
	    		// add to the bucket (if size allows)
	    		uint32 curr_size = buckets->per_thread_bucket_sizes[bucket_index][tid];
	    		if(curr_size + 1 <= params->bucket_size) {
	    			bool store_pos = true;
	    			// don't store if near-by sequence present, need to check the last value only
	    			if(buckets->per_thread_bucket_sizes[bucket_index][tid] > 0) {
	    				seq_t L = pos > params->bucket_entry_coverage ? pos - params->bucket_entry_coverage : 0;
	    				seq_t epos = bucket[buckets->per_thread_bucket_sizes[bucket_index][tid]-1];
	    				if(epos >= L) {
	    					store_pos = false;
	    				}
	    			}
	    			if(store_pos) {
	    				bucket[curr_size] = pos;
	    				buckets->per_thread_bucket_sizes[bucket_index][tid]++;
	    				n_bucket_entries++;
	    			} else {
	    				n_filtered++;
	    			}
	    		} else {
	    			n_dropped++;
	    		}
	    	}
	    }
	}
	end_time = omp_get_wtime();
	printf("Populated all the buckets. Time : %.2f sec\n", end_time - start_time);

	// 4. collect per thread results
	printf("Collecting thread buckets... \n");
	double start_coll_sort = omp_get_wtime();
	for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
		buckets_t* buckets = &ref.hash_tables[t];
		for(uint32 b = 0; b < buckets->next_free_bucket_index; b++) {
			VectorSeqPos& bucket = buckets->buckets_data_vectors[b];
			for(uint32 i = 0; i < params->n_threads; i++) {
				VectorSeqPos& thread_bucket = buckets->per_thread_buckets_data_vectors[b][i];
				for(uint32 j = 0; j < buckets->per_thread_bucket_sizes[b][i]; j++) {
					bucket.push_back(thread_bucket[j]);
					buckets->bucket_sizes[b]++;
				}
				thread_bucket = VectorSeqPos();
			}

		}
	}
	end_time = omp_get_wtime();
	printf("Collected all the buckets. Time : %.2f sec\n", end_time - start_coll_sort);

	// 5. sort each bucket!
	printf("Sorting buckets... \n");
	double start_time_sort = omp_get_wtime();
	#pragma omp parallel for
	for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
		buckets_t* buckets = &ref.hash_tables[t];
		for(uint32 b = 0; b < buckets->next_free_bucket_index; b++) {
			VectorSeqPos& bucket = buckets->buckets_data_vectors[b];
			std::sort(bucket.begin(), bucket.begin() + buckets->bucket_sizes[b]);
		}
	}
	end_time = omp_get_wtime();
	printf("Total sort time : %.2f sec\n", end_time - start_time_sort);

	printf("Total number of valid reference windows: %u \n", n_valid_windows);
	printf("Total number of valid reference windows with valid hashes: %u \n", n_valid_hashes);
	printf("Total number of window bucket entries: %llu \n", n_bucket_entries);
	printf("Total number of window bucket entries filtered: %llu \n", n_filtered);
	printf("Total number of window bucket entries dropped: %llu \n", n_dropped);
	end_time = omp_get_wtime();
	printf("Total hashing time: %.2f sec\n", end_time - start_time);
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
	//compute_kmer_counts(ref.seq.c_str(), ref.seq.size(), params, ref.kmer_hist);
	//store_kmer_hist(fastaFname, ref.kmer_hist);
	//ref.kmer_hist = MapKmerCounts(); // free memory
	load_freq_kmers(fastaFname, ref.high_freq_kmer_trie, params->max_count);
	mark_freq_kmers(ref, params);
	double end_time = omp_get_wtime();
	printf("Total kmer pre-processing time: %.2f sec\n", end_time - start_time);

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
