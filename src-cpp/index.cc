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
			ref.windows_by_pos.insert(std::pair<seq_t, ref_win_t>(i, window));
		} else {
			// skip until at least 1 potential valid kmer
			//i += params->k;
		}
	}
}

// --- Simhash / Minhash ---

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
	
	// 2. compute the frequency of each kmer
	//	  and collect high-frequency kmers
	t = clock();
	if(params->kmer_type != SPARSE) {
		compute_kmer_counts(ref.seq.c_str(), ref.seq.size(), params, ref.kmer_hist);
		find_high_freq_kmers(ref.kmer_hist, ref.high_freq_kmer_hist, params);
		ref.kmer_hist = MapKmerCounts(); // free memory
		printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	
	// 3. compute valid reference windows
	t = clock();
	generate_ref_windows(ref, params);
	printf("Total number of valid windows: %zu\n", ref.windows_by_pos.size());
	
	// 4. hash each window

	if (params->alg == MINH) {
		// initialize the hash tables
		//ref.minhash_maps_by_h.resize(params->h);
		ref.hash_buckets.resize(params->h/params->band_size);
		printf("Total number of buckets: %u \n", params->h/params->band_size);
	}

	// collect the valid reference positions for parallel iteration
	VectorSeqPos valid_positions;
	for(MapPos2Window::iterator it = ref.windows_by_pos.begin(); it != ref.windows_by_pos.end(); it++) {
		valid_positions.push_back(it->first);
	}

	t = clock();
	#pragma omp parallel for
	for(seq_t i = 0; i < valid_positions.size(); i++) {
		ref_win_t* w = &ref.windows_by_pos[valid_positions[i]];
		if(params->alg == SIMH) {
			w->simhash = simhash(ref.seq.c_str(), w->pos, params->ref_window_size,
					ref.high_freq_kmer_hist, MapKmerCounts(), params, 1);
		} else if (params->alg == MINH) {
			w->minhashes.resize(params->h);
			w->simhash = minhash(ref.seq.c_str(), w->pos, params->ref_window_size,
					ref.high_freq_kmer_hist, MapKmerCounts(), params, 1, w->minhashes);

			for(uint32 band = 0; band < params->h/params->band_size; band++) {
				std::string band_entries;
				for(uint32 v = 0; v < params->band_size; v++) {
					band_entries += std::to_string(w->minhashes[params->band_size*band + v]);
					band_entries += std::string(".");
				}
				minhash_t hash = CityHash32(band_entries.c_str(), band_entries.size());

				#pragma omp critical
				{
					std::map<minhash_t, VectorSeqPos>::iterator v;
					if((v = ref.hash_buckets[band].find(hash)) != ref.hash_buckets[band].end()) {
						v->second.push_back(w->pos);
					} else {
						ref.hash_buckets[band].insert(std::pair<minhash_t, VectorSeqPos>(hash, VectorSeqPos()));
						ref.hash_buckets[band][hash].push_back(w->pos);
					}
				}
			}

			// insert window into the minhash maps based on its minhash value for each hash function
//			for(uint32 h = 0; h < params->h; h++) {
//				minhash_t minh = w->minhashes[h];
//
//				#pragma omp critical
//				{
//					std::map<minhash_t, VectorWindowPtr>::iterator v;
//					if((v = ref.minhash_maps_by_h[h].find(minh)) != ref.minhash_maps_by_h[h].end()) {
//						v->second.push_back(w);
//					} else {
//						ref.minhash_maps_by_h[h].insert(std::pair<minhash_t, VectorWindowPtr>(minh, VectorWindowPtr()));
//						ref.minhash_maps_by_h[h][minh].push_back(w);
//					}
//				}
//			}
			w->minhashes = VectorMinHash();
		}
	}
	printf("Total hashing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	if (params->alg == SIMH) {
		// 6. sort
		t = clock();
		sort_windows_hash(ref);
		printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		printf("Total number of distinct window hashes = %llu \n", (uint64) get_num_distinct(ref));
	}
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
			read_t r = reads.reads[i];
			compute_kmer_counts(r.seq.c_str(), r.len, params, reads.kmer_hist);
		}
		find_low_freq_kmers(reads.kmer_hist, reads.low_freq_kmer_hist, params);
		reads.kmer_hist = MapKmerCounts(); // free memory
		printf("Total kmer histogram generation time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}

	// 3. compute the fingerprints of each read
	t = clock();

	#pragma omp parallel for
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(params->alg == SIMH) {
			r->simhash = simhash(r->seq.c_str(), 0, r->len, ref.high_freq_kmer_hist, reads.low_freq_kmer_hist, params, 0);
		} else if (params->alg == MINH) {
			r->minhashes.resize(params->h);
			r->simhash = minhash(r->seq.c_str(), 0, r->len, ref.high_freq_kmer_hist, reads.low_freq_kmer_hist, params, 0, r->minhashes);
		}
	}
	printf("Total hashing time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

	// free memory
	ref.high_freq_kmer_hist = MapKmerCounts();
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
void index_ref_table_i(ref_t& ref, const index_params_t* params, const seq_t i) {
	// hash each window using sampling ids i
	clock_t t = clock();
	for(seq_t j = 0; j < ref.windows_by_pos.size(); j++) {
		ref.windows_by_pos[j].simhash = sampling_hash(ref.seq.c_str(), ref.windows_by_pos[j].pos, i, params);
	}
	printf("Total hash table %llu computation time: %.2f sec\n", (uint64) i, (float)(clock() - t) / CLOCKS_PER_SEC);
	
	// sort the hashes
	t = clock();
	sort_windows_hash(ref);
	printf("Total sorting time: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	printf("Total number of distinct window hashes = %llu \n", (uint64) get_num_distinct(ref));
}

void index_reads_table_i(reads_t& reads, const index_params_t* params, const seq_t i) {
	// hash each read using sampling ids i
	clock_t t = clock();
	for(uint32 j = 0; j < reads.reads.size(); j++) {
		reads.reads[j].simhash = sampling_hash(reads.reads[j].seq.c_str(), 0, i, params);
	}
	printf("Total read hash table %llu computation time: %.2f sec\n", (uint64) i, (float)(clock() - t) / CLOCKS_PER_SEC);
}
