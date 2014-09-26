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

// --- Minhash ---

// index reference:
// - load fasta
// - generate valid windows (sliding window)
// - compute the hash of each window
// - sort
void index_ref_lsh(const char* fastaFname, index_params_t* params, ref_t& ref) {
	printf("**** SRX Reference Indexing ****\n");
		
	// 1. load the reference
	clock_t t = clock();
	fasta2ref(fastaFname, ref);
	printf("Total ref loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// 2. load the frequency of each kmer and collect high-frequency kmers
	t = clock();
	if(params->max_count != 0) {
		//compute_kmer_counts(ref.seq.c_str(), ref.seq.size(), params, ref.kmer_hist);
		//store_kmer_hist(fastaFname, ref.kmer_hist);
		//load_kmer_hist(ref.seq.c_str(), ref.kmer_hist, params->max_count);
		//find_high_freq_kmers(ref.kmer_hist, ref.high_freq_kmer_hist, params);

		load_freq_kmers(fastaFname, ref.high_freq_kmer_trie, params->max_count);
		ref.kmer_hist = MapKmerCounts(); // free memory
		printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	
	// 3. hash each valid window

	// initialize the hash tables
	ref.hash_tables.resize(params->n_tables);
	for(uint32 t = 0; t < params->n_tables; t++) {
		buckets_t* buckets = &ref.hash_tables[t];
		buckets->n_buckets = pow(2, params->n_buckets_pow2);
		buckets->next_free_bucket_index = 0;
		//buckets->bucket_sizes.resize(buckets->n_buckets);
		buckets->bucket_indices.resize(buckets->n_buckets, buckets->n_buckets);
		buckets->buckets_data_vectors.resize(buckets->n_buckets);
	}

	uint32 n_valid_windows = 0;
	uint32 n_bucket_entries = 0;
	uint32 n_filtered = 0;

	t = clock();
	#pragma omp parallel for
	for(seq_t pos = 0; pos < ref.seq.size() - params->ref_window_size + 1; pos++) { // for each window of the genome
		if(!is_inform_ref_window(&ref.seq.c_str()[pos], params->ref_window_size)) {
			continue; // discard windows with low information content
		}

		#pragma omp atomic
		n_valid_windows++;

		VectorMinHash minhashes(params->h); // TODO: each thread should index into its pre-allocated buffer
		// get the min-hash signature for the window
		minhash(ref.seq.c_str(), pos, params->ref_window_size, ref.high_freq_kmer_trie,  marisa::Trie(), params, 1, minhashes);

		for(uint32 t = 0; t < params->n_tables; t++) { // for each hash table
			VectorMinHash sketch_proj(params->sketch_proj_len);
			for(uint32 p = 0; p < params->sketch_proj_len; p++) {
				sketch_proj[p] = minhashes[params->sketch_proj_indices[t*params->sketch_proj_len + p]];
			}
			minhash_t bucket_hash = params->sketch_proj_hash_func.apply_vector(sketch_proj);

			#pragma omp critical
			{
				buckets_t* buckets = &ref.hash_tables[t];
				uint32 bucket_index = buckets->bucket_indices[bucket_hash];
				if(bucket_index == buckets->n_buckets) {
					// this is the first entry in the bucket
					buckets->bucket_indices[bucket_hash] = buckets->next_free_bucket_index;
					buckets->next_free_bucket_index++;

					VectorSeqPos& bucket = buckets->buckets_data_vectors[buckets->bucket_indices[bucket_hash]];
					//bucket.reserve(params->bucket_size);
					bucket.push_back(pos); // store the window position in the bucket
					//buckets->bucket_sizes[bucket_hash]++;
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
						//buckets->bucket_sizes[bucket_hash]++;
					}
				}
			}
		}
	}
	printf("Total number of valid reference windows: %u \n", n_valid_windows);
	printf("Total number of window bucket entries: %u \n", n_bucket_entries);
	printf("Total number of window bucket entries filtered: %u \n", n_filtered);
	printf("Total hashing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
}

void index_reads_lsh(const char* readsFname, ref_t& ref, index_params_t* params, reads_t& reads) {
	printf("**** SRX Read Indexing ****\n");

	// 1. load the reads (TODO: batch mode)
	clock_t t = clock();
	fastq2reads(readsFname, reads);
	printf("Total read loading time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// 2. compute the frequency of each kmer
	t = clock();
	if(params->kmer_type != SPARSE) {
		for(uint32 i = 0; i < reads.reads.size(); i++) { // add the contribution of each read
			//read_t r = reads.reads[i];
			//compute_kmer_counts(r.seq.c_str(), r.len, params, reads.kmer_hist);
		}
		//find_low_freq_kmers(reads.kmer_hist, reads.low_freq_kmer_hist, params);
		reads.kmer_hist = MapKmerCounts(); // free memory
		printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	// 3. compute the fingerprints of each read
	t = clock();

	#pragma omp parallel for
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		r->minhashes.resize(params->h);
		minhash(r->seq.c_str(), 0, r->len, ref.high_freq_kmer_trie,  marisa::Trie(), params, 0, r->minhashes);
	}
	printf("Total hashing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// free memory
//	ref.high_freq_kmer_trie =  marisa::Trie();
	reads.low_freq_kmer_hist = MapKmerCounts();

	// 4. sort the reads by their hash
	/*sort_reads_hash(reads);

	// 5. split reads into "clusters" based on their simhash
	t = clock();
	clusters_t* clusters;
	cluster_sorted_reads(reads, &clusters);
	printf("Total number of clusters = %d \n", clusters->num_clusters);
	printf("Total clustering time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	//int count = 0;
	//for(int i = 0; i < clusters->num_clusters; i++) {
		//if((count < 20) && (clusters->clusters[i].size > 1)) {
			//printf("cluster = %d, simhash = %llx, size = %d \n", i, clusters->clusters[i].simhash, clusters->clusters[i].size);
			//count++;
			//for(int j = 0; j < clusters->clusters[i].size; j++) {
				//print_read(clusters->clusters[i].reads[j]);
			//}
		//}
	//}

	// 6. collapse the clusters that are close to each other in Hamming distance
	t = clock();
	int num_collapsed = collapse_clusters(clusters, params);
	printf("Collapsed %d clusters \n", num_collapsed);
	printf("Total number of clusters remaining = %d \n", clusters->num_clusters - num_collapsed);
	printf("Total clustering collapse time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	*/
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
