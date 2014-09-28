#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include <time.h>
#include <vector>
#include <string.h>

#include "index.h"
#include "io.h"
#include "hash.h"
#include "cluster.h"

// generate all valid/informative reference windows
void generate_ref_windows(ref_t& ref, index_params_t* params) {
	seq_t max_num_windows = ref.seq.size() - params->ref_window_size + 1;
	seq_t i;
	for(i = 0; i < max_num_windows; i++) {
		if(is_inform_ref_window(&ref.seq.c_str()[i], params->ref_window_size)) {
			ref_win_t window;
			window.pos = i;
			//ref.windows_by_pos.insert(std::pair<seq_t, ref_win_t>(i, window));
		}
	}
}

void mark_freq_kmers(ref_t& ref, const index_params_t* params) {
	clock_t t = clock();
	ref.ignore_kmer_bitmask.resize(ref.len - params->k);
	marisa::Agent agent;
	for(seq_t i = 0; i < ref.len - params->k; i++) {
		for (uint32 k = 0; k < params->k; k++) {
			if(ref.seq.c_str()[k] == BASE_IGNORE) {
				ref.ignore_kmer_bitmask[i] = 1; // contains ambiguous bases
				break;
			}
		}
		agent.set_query(&ref.seq.c_str()[i], params->k);
		if(ref.high_freq_kmer_trie.lookup(agent)) {
			ref.ignore_kmer_bitmask[i] = 1;
		}
	}
	printf("Done marking frequent kmers time: %.2f sec \n", (float) (clock() - t)/CLOCKS_PER_SEC);
}

void mark_windows_to_discard(ref_t& ref, const index_params_t* params) {
	ref.ignore_window_bitmask.resize(ref.len- params->ref_window_size + 1);
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

	// 2. mark uninformative windows / load
	printf("Marking uninformative windows... \n");
	t = clock();
	mark_windows_to_discard(ref, params);
	printf("Done marking windows; time: %.2f sec \n", (float) (clock() - t)/CLOCKS_PER_SEC);
	
	// 3. load the frequency of each kmer and collect high-frequency kmers
	t = clock();
	//compute_kmer_counts(ref.seq.c_str(), ref.seq.size(), params, ref.kmer_hist);
	//store_kmer_hist(fastaFname, ref.kmer_hist);
	//ref.kmer_hist = MapKmerCounts(); // free memory
	load_freq_kmers(fastaFname, ref.high_freq_kmer_trie, params->max_count);
	mark_freq_kmers(ref, params);
	printf("Total kmer pre-processing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 4. hash each valid window
	ref.hash_tables.resize(params->n_tables); // initialize the hash tables
	for(uint32 t = 0; t < params->n_tables; t++) {
		buckets_t* buckets = &ref.hash_tables[t];
		buckets->n_buckets = pow(2, params->n_buckets_pow2);
		buckets->next_free_bucket_index = 0;
		buckets->bucket_indices.resize(buckets->n_buckets, buckets->n_buckets);
		buckets->buckets_data_vectors.resize(buckets->n_buckets);
	}

	// per-thread storage
	std::vector<VectorMinHash> minhash_thread_vectors(params->n_threads);
	for(uint32 i = 0; i < params->n_threads; i++) {
		minhash_thread_vectors[i].resize(params->h);
		std::fill(minhash_thread_vectors[i].begin(), minhash_thread_vectors[i].end(), UINT_MAX);
	}
	std::vector<CyclicHash*> hasher_thread_vectors(params->n_threads);
	for(uint32 i = 0; i < params->n_threads; i++) {
		hasher_thread_vectors[i] = new CyclicHash(params->k, 32);
	}

	printf("Hashing reference windows... \n");
	uint32 n_valid_windows = 0;
	uint32 n_valid_hashes = 0;
	uint32 n_bucket_entries = 0;
	uint32 n_filtered = 0;
	t = clock();
	omp_set_num_threads(params->n_threads); // split the windows across the threads
	#pragma omp parallel reduction(+:n_valid_windows, n_valid_hashes, n_bucket_entries, n_filtered)
	{
	    int tid = omp_get_thread_num();
	    int n_threads = omp_get_num_threads();
	    seq_t chunk_start = tid*(ref.len - params->ref_window_size + 1) / n_threads;
	    seq_t chunk_end = (tid + 1)*(ref.len - params->ref_window_size + 1) / n_threads;

	    for (seq_t pos = chunk_start; pos != chunk_end; pos++) { // for each window of the thread's chunk
	    	/*if(ref.ignore_window_bitmask[pos]) { // discard windows with low information content
				continue;
			}*/
	    	n_valid_windows++;
	    	if((n_valid_windows) % 1000000 == 0) {
	    		printf("Thread %d processed %u valid windows \n", tid, n_valid_windows);
	    	}

	    	// get the min-hash signature for the window
	    	VectorMinHash& minhashes = minhash_thread_vectors[0]; // each thread indexes into its pre-allocated buffer
	    	bool valid_hash;
	    	if(pos == chunk_start) {
	    		valid_hash = minhash_rolling_init(ref, pos, params->ref_window_size,
							ref.ignore_kmer_bitmask, params,
							hasher_thread_vectors[0],
							minhashes);
	    	} else {
	    		valid_hash = minhash_rolling(ref, pos, params->ref_window_size,
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
				#pragma omp critical
	    		{
	    			buckets_t* buckets = &ref.hash_tables[t];
	    			uint32 bucket_index = buckets->bucket_indices[bucket_hash];
	    			if(bucket_index == buckets->n_buckets) {
	    				// this is the first entry in the bucket
	    				buckets->bucket_indices[bucket_hash] = buckets->next_free_bucket_index;
	    				buckets->next_free_bucket_index++;
	    				VectorSeqPos& bucket = buckets->buckets_data_vectors[buckets->bucket_indices[bucket_hash]];
	    				bucket.push_back(pos); // store the window position in the bucket
	    			} else {
	    				VectorSeqPos& bucket = buckets->buckets_data_vectors[bucket_index];
	    				// add to the existing hash bucket
	    				if(bucket.size() + 1 < params->bucket_size) {
	    					// don't store if near-by sequence present
	    					bool store_pos = true;
	    					for(uint32 e = 0; e < bucket.size(); e++) {
	    						seq_t epos = bucket[e];
	    						seq_t H = pos + params->bucket_entry_coverage;
	    						seq_t L = pos > params->bucket_entry_coverage ? pos - params->bucket_entry_coverage : 0;
	    						if((epos <= H) && (epos >= L)) {
	    							store_pos = false;
	    							break;
	    						}
	    					}
	    					if(store_pos) {
	    						bucket.push_back(pos);
	    						n_bucket_entries++;
	    					} else {
	    						n_filtered++;
	    					}
	    				}
	    			}
	    		}
	    	}
	    }
	}

	printf("Total number of valid reference windows: %u \n", n_valid_windows);
	printf("Total number of valid reference windows with valid hashes: %u \n", n_valid_hashes);
	printf("Total number of window bucket entries: %u \n", n_bucket_entries);
	printf("Total number of window bucket entries filtered: %u \n", n_filtered);
	printf("Total hashing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
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
	t = clock();
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
	printf("Total read hashing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
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
