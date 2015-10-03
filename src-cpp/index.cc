#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <omp.h>
#include "types.h"
#include "index.h"
#include "lsh.h"
#include "io.h"
#include "hash.h"
#include <fstream>


// compute and store the frequency of each kmer in the given sequence (up to length 16)
void compute_kmer_hist16(const char* seq, const seq_t seq_len, const index_params_t* params, MapKmerCounts& hist) {
	for(seq_t j = 0; j <= (seq_len - params->k); j++) {
		uint32_t kmer;
		if(pack_32(&seq[j], params->k, &kmer) < 0) {
			continue;
		}
		hist[kmer]++;
	}
}


void ref_kmer_fingerprint_stats(const char* fastaFname, index_params_t* params, ref_t& ref) {
	printf("Loading FASTA file %s... \n", fastaFname);
	clock_t t = clock();
	fasta2ref(fastaFname, ref);
	printf("Reference loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	std::vector<std::string> freq_fingerprints(ref.len - params->ref_window_size + 1);
	uint32 n_kmers = (params->ref_window_size - params->k2 + 1);
	
	uint32 unique = 0; 
	uint32 repk = 0;
	#pragma omp parallel for reduction(+:unique, repk)
	for(seq_t pos = 0; pos < ref.len - params->ref_window_size + 1; pos++) { // for each window of the genome
		// generate the kmers of this window (pairs of kmers + pos)
		std::vector<std::pair<minhash_t, seq_t>> kmers(n_kmers);
		for(uint32 j = 0; j < n_kmers; j++) {
			kmers[j] = std::make_pair(CityHash32(&ref.seq.c_str()[pos + j], params->k2), j);
		}
		// sort the kmers
		std::sort(kmers.begin(), kmers.end());

		// unique = -1
		// repeat = (list of positions)
		std::string window;
		//window += std::to_string(pos);
		//window += ": ";
		std::vector<uint16_t> counts;
		uint16_t counter = 0;
		for(uint32 j = 0; j < n_kmers; j++) {
			if(j == 0) {
				counter = 1;
			} else {
				if(kmers[j].first == kmers[j-1].first) {
					counter++;
					if(j == n_kmers -1) {
						counts.push_back(counter);
						//for(int k = 0; k < counter; k++) {
						//	window += std::to_string(kmers[j-k].second);
						//	window += ",";
						//}
						//window += ";";
					}
				} else {
					if(counter > 1) {
						counts.push_back(counter);
						//if(counter > 2) break;	
						//for(int k = 0; k < counter; k++) {
						//	window += kmers[j-1-k].second;
						//	window += ",";
						//}
						//window += ";";
					}
					counter = 1;
				}
			}
			//if(counts.size() > 0) break; // TEMP
		}
		if(counts.size() > 0) {
			//if(iupacChar[(int)ref.seq[pos]] == 'N') continue;
			//for(uint32 x = 0; x < params->ref_window_size; x++) {
                        //         printf("%c", iupacChar[(int)ref.seq[pos + x]]);
                        //}
			
			std::sort(counts.begin(), counts.end());
			for(int k = 0; k < counts.size(); k++) {
				//if(counts[k] > 2) {
				//	repk += 1;
				//	break;
				//}
				window += std::to_string(counts[k]);
				window += " ";
			} 
			//std::cout << "\n";
			//std::cout<<window;
		} else {
			unique += 1;
			//window += "-1";
		}
		freq_fingerprints[pos] = window;
		//if(pos % 10000000 == 0) {
		//	std::cout << "Processed " << pos << " windows\n";
		//}
	}

	
	std::cout << "Computed fingerprints, unique " << unique << " repeats > 2 " << repk << "\n";

	std::string stats_fname(fastaFname);
        stats_fname += std::string("__fingerprints");
        std::ofstream stats_file;
        stats_file.open(stats_fname.c_str(), std::ios::out);
        if (!stats_file.is_open()) {
                printf("store_ref_idx: Cannot open the KMER FINGERPRINT file %s!\n", stats_fname.c_str());
                exit(1);
        }
	std::cout << "Total number of fingerprints: " << freq_fingerprints.size() << "\n";
	std::sort(freq_fingerprints.begin(), freq_fingerprints.end());
	int counter = 0;
	int repeat = 0;
	for(uint32 j = 0; j < freq_fingerprints.size(); j++) {
		if(j == 0) {
			if(freq_fingerprints[j] != "-1") {
				counter = 1;
			}
		} else {
			if(freq_fingerprints[j] != "-1" && freq_fingerprints[j] == freq_fingerprints[j-1]) {
				counter++;
				if(j == freq_fingerprints.size()-1) {
					repeat += counter;
					stats_file << freq_fingerprints[j] << "," << counter << "\n";
				}
			} else {
				if(counter >= 1) {
					repeat += counter;
					stats_file << freq_fingerprints[j-1] << "," << counter << "\n";
				}
				if(freq_fingerprints[j] != "-1") {
					counter = 1;
				} else {
					counter = 0;
				}
			}
		}
	}
	std::cout << "Total number of unique fingerprints: " << freq_fingerprints.size() - repeat << "\n";

	//for(uint32 i = 0; i < freq_fingerprints.size(); i++) {
	//	stats_file << freq_fingerprints[i]  << "\n";
	//}
	stats_file.close();
}

void mark_freq_kmers(ref_t& ref, const index_params_t* params) {
	clock_t t = clock();
	ref.ignore_kmer_bitmask.resize(ref.len - params->k);
	#pragma omp parallel for
	for(seq_t i = 0; i < ref.len - params->k + 1; i++) { // for each window of the genome
#if USE_MARISA
		marisa::Agent agent;
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
#else
		uint32_t packed_kmer;
		if(pack_32(&ref.seq.c_str()[i], params->k, &packed_kmer) < 0) {
			ref.ignore_kmer_bitmask[i] = 1; // contains ambiguous bases
			continue;
		}
		if(ref.high_freq_kmer_bitmap[packed_kmer]) {
			ref.ignore_kmer_bitmask[i] = 1;
		}
#endif
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
// - compute the fingerprint of each window
// - bucket into multiple hash tables using different projections
// - sort buckets

#define REPORT_WINDOW_PROC_GRANULARITY 2000000

void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref) {
	// 1. load the reference
	printf("Loading FASTA file %s... \n", fastaFname);
	clock_t t = clock();
	fasta2ref(fastaFname, ref);
	printf("Reference loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 2. load the frequency of each kmer and collect high-frequency kmers
	printf("Loading frequent kmers... \n");
	double start_time = omp_get_wtime();
	//compute_kmer_hist16(ref.seq.c_str(), ref.seq.size(), params, ref.kmer_hist);
	//store_kmer_hist(fastaFname, ref.kmer_hist);
	//ref.kmer_hist = MapKmerCounts(); // free memory
	load_freq_kmers(fastaFname, ref.high_freq_kmer_bitmap, ref.high_freq_kmer_trie, params->max_count);
	mark_freq_kmers(ref, params);

	printf("Loading valid windows mask... \n");
	if(!load_valid_window_mask(fastaFname, ref, params)) {
		mark_windows_to_discard(ref, params);
		store_valid_window_mask(fastaFname, ref, params);
	}
	printf("Total window/kmer pre-processing time: %.2f sec\n", omp_get_wtime() - start_time);

	// initialize the hash tables
	ref.mutable_index.per_table_buckets.resize(params->n_tables);
	for(uint32 t = 0; t < params->n_tables; t++) {
		buckets_t* buckets = &ref.mutable_index.per_table_buckets[t];
		// initialize global buckets
		buckets->buckets_data_vectors.resize(params->n_buckets);
		// initialize per thread buckets
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

	// 3. hash each valid window
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
	    seq_t chunk_start = ((ref.len - params->ref_window_size + 1) / n_threads)*tid;
	    seq_t chunk_end = ((ref.len - params->ref_window_size + 1) / n_threads)*(tid + 1);
	    printf("Thread %d range: %u %u \n", tid, chunk_start, chunk_end);

	    int sync_point = 1;
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
	    	if(ref.ignore_window_bitmask[pos]) {
	    		init_minhash = true;
	    		continue;
			}
	    	n_valid_windows++;

	    	// get the min-hash signature for the window
	    	VectorMinHash& minhashes = minhash_thread_vectors[tid]; // each thread indexes into its pre-allocated buffer
	    	minhash_matrix_t& rolling_minhash_matrix = minhash_matrices[tid];
	    	bool valid_hash;
	    	if(init_minhash == true) {
	    		valid_hash = minhash_rolling_init(ref.seq.c_str(), pos, params->ref_window_size,
	    					rolling_minhash_matrix, ref.ignore_kmer_bitmask, params,
							minhashes);
	    		init_minhash = false;
	    	} else {
	    		valid_hash = minhash_rolling(ref.seq.c_str(), pos, params->ref_window_size,
	    					rolling_minhash_matrix, ref.ignore_kmer_bitmask, params,
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
	    		buckets_t* buckets = &ref.mutable_index.per_table_buckets[t];
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
					if((epos->pos + epos->len) == pos) {
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
		buckets_t* buckets = &ref.mutable_index.per_table_buckets[t];
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
		buckets_t* buckets = &ref.mutable_index.per_table_buckets[t];
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

void load_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref) {
	printf("Loading FASTA file %s... \n", fastaFname);
	clock_t t = clock();
	fasta2ref(fastaFname, ref);
	printf("Reference loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	printf("Loading reference index for reference file %s... \n", fastaFname);
	t = clock();
	load_ref_idx(fastaFname, ref, params);
	load_freq_kmers(fastaFname, ref.high_freq_kmer_bitmap, ref.high_freq_kmer_trie, params->max_count);
	printf("Reference index loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// load/precompute k2 voting kmers
	t = clock();
	printf("Precomputing/loading k2 kmer hashes... \n");
	if(params->precomp_k2 && !load_kmer2_hashes(fastaFname, ref, params)) {
		omp_set_num_threads(params->n_threads);
		compute_store_kmer2_hashes(fastaFname, ref, params);
	}
	printf("Reference k2 hash loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}

void store_index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref) {
	printf("Storing the reference index for reference file %s... \n", fastaFname);
	clock_t t = clock();
	store_ref_idx(fastaFname, ref, params);
	printf("Reference index storing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}

void index_reads_lsh(const char* readsFname, ref_t& ref, index_params_t* params, reads_t& reads) {
	// load the reads (TODO: batch mode)
	printf("Loading FATQ reads...\n");
	clock_t t = clock();
	fastq2reads(readsFname, reads);
	printf("Total read loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// compute the minhash fingerprint of each read
	printf("Hashing reads...\n");
	double start_time = omp_get_wtime();
	omp_set_num_threads(params->n_threads);
	#pragma omp parallel for
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		r->minhashes.resize(params->h);
		r->valid_minhash = minhash(r->seq.c_str(), r->len,
				ref.high_freq_kmer_bitmap,
				ref.high_freq_kmer_trie,
				MarisaTrie(), params,
				r->minhashes);

		r->minhashes_rc.resize(params->h);
		r->valid_minhash_rc = minhash(r->rc.c_str(), r->len,
				ref.high_freq_kmer_bitmap,
				ref.high_freq_kmer_trie,
				MarisaTrie(), params,
				r->minhashes_rc);
	}
	double end_time = omp_get_wtime();
	printf("Total read hashing time: %.2f sec\n", end_time - start_time);
}
