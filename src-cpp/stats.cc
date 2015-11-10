#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <time.h>
#include <omp.h>
#include <fstream>
#include <limits.h>
#include "types.h"
#include "io.h"
#include "hash.h"

// K-Mer stats //

void load_freq_kmers(const char* refFname, std::set<uint64>& freq_kmers, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".kmer_hist");
	fname += std::to_string(params->k);
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("load_freq_kmers: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}
	uint32 filtered = 0;
	uint64 kmer;
	int count;
	while(true) {
		file.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
		if(!file) {
			break;
		}
		file.read(reinterpret_cast<char*>(&count), sizeof(count));
		if(count >= params->max_count) {
			freq_kmers.insert(kmer);
			filtered++;
		}
	}
	file.close();
	std::set<uint64>::iterator it;
	printf("Filtered %u kmers\n", filtered);
}

void load_freq_kmers(const char* refFname, VectorBool& freq_kmers_bitmap, MarisaTrie& freq_trie, const uint32 max_count_threshold) {
	std::string fname(refFname);
	fname += std::string(".kmer_hist");

	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);

	if (!file.is_open()) {
		printf("load_kmer_hist: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}
	freq_kmers_bitmap.resize(UINT_MAX + 1L);
	uint32 kmer, count;
	uint32 filtered = 0;
	uint32 tot_filtered = 0;
	int map_size;
	file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
#if(USE_MARISA)
	marisa::Keyset keys;
#endif
	while(map_size >= 0) {
		file.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
		file.read(reinterpret_cast<char*>(&count), sizeof(count));
		if(count >= max_count_threshold) {
			tot_filtered++;
			//if (!kmer_has_zero(kmer)) {
			freq_kmers_bitmap[kmer] = true;
			filtered++;
			//}
#if(USE_MARISA)
			unsigned char* seq = (unsigned char*) malloc(17*sizeof(char));
			unpack_32(kmer, seq, 16);
			seq[16] = '\0';
			keys.push_back((const char*) seq);
			free(seq);
#endif
		}
		map_size--;
	}
#if(USE_MARISA)
	freq_trie.build(keys, 0);
#endif
	file.close();
	printf("Filtered %u kmers, tot %u \n", filtered, tot_filtered);
}

// compute and store the frequency of each kmer in the given sequence (up to length 32)
void compute_and_store_kmer_hist32(const char* refFname, const char* seq, const seq_t seq_len, const index_params_t* params) {
	std::vector<uint64> kmers(seq_len - params->k + 1);
	#pragma omp parallel
	for(seq_t j = 0; j < seq_len - params->k + 1; j++) {
			pack_64(&seq[j], params->k, &(kmers[j]));
	}
	std::sort(kmers.begin(), kmers.end());
	std::string fname(refFname);
	fname += std::string(".kmer_hist");
	fname += std::to_string(params->k);
	std::ofstream file;
    file.open(fname.c_str(), std::ios::out | std::ios::binary);
    if (!file.is_open()) {
		printf("compute_and_store_kmer_hist32: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}

    int counter = 0;
    seq_t unique = 0;
	for(seq_t j = 0; j < kmers.size(); j++) {
		if(j == 0) {
			if(kmers[j] != 0) {
				counter = 1;
			}
		} else {
			if(kmers[j] != 0 && kmers[j] == kmers[j-1]) {
				counter++;
				if(j == kmers.size()-1) {
					file.write(reinterpret_cast<char*>(&kmers[j]), sizeof(kmers[j]));
					file.write(reinterpret_cast<char*>(&counter), sizeof(counter));
				}
			} else {
				if(counter >= 1) {
					if(counter == 1) {
						unique++;
					} else {
						file.write(reinterpret_cast<char*>(&kmers[j-1]), sizeof(kmers[j-1]));
						file.write(reinterpret_cast<char*>(&counter), sizeof(counter));
					}
				}
				if(kmers[j] != 0) {
					counter = 1;
				} else {
					counter = 0;
				}
			}
		}
	}
	file.close();
    std::cout << "Total number of unique kmers " << unique << " out of " << kmers.size() << "\n";
}

// compute and store the frequency of each kmer in the given sequence (up to length 16)
void compute_and_store_kmer_hist16(const char* refFname, const char* seq, const seq_t seq_len, const index_params_t* params) {
	std::map<uint32, seq_t> hist;
	for(seq_t j = 0; j < seq_len - params->k + 1; j++) {
		uint32_t kmer;
		if(pack_32(&seq[j], params->k, &kmer) < 0) {
			continue;
		}
		hist[kmer]++;
	}
	std::string fname(refFname);
	fname += std::string(".kmer_hist");
	fname += std::to_string(params->k);
	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::binary);
	if (!file.is_open()) {
		printf("compute_and_store_kmer_hist16: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}
	uint32 map_size = hist.size();
	file.write(reinterpret_cast<char*>(&map_size), sizeof(map_size));
	for (std::map<uint32, seq_t>::const_iterator it = hist.begin(); it != hist.end(); ++it) {
		file.write(reinterpret_cast<const char*>(&(it->first)), sizeof(it->first));
		file.write(reinterpret_cast<const char*>(&(it->second)), sizeof(it->second));
	}
	file.close();
}

#define NUM_HIST_BUCKETS 1000000
void store_kmer_hist_stat(const char* refFname, const MapKmerCounts& hist) {
	std::string fname(refFname);
	fname += std::string(".kmer_hist_stats");
	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::app);
	if (!file.is_open()) {
		printf("store_kmer_hist: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}
	VectorU32 freq_buckets(NUM_HIST_BUCKETS);
	for (MapKmerCounts::const_iterator it = hist.begin(); it != hist.end(); ++it) {
		//uint32 kmer = it->first;
		seq_t count = it->second;
		if(count >= NUM_HIST_BUCKETS) {
			freq_buckets[NUM_HIST_BUCKETS-1]++;
		} else {
			freq_buckets[count]++;
		}
	}
	for(uint32 i = 0; i < NUM_HIST_BUCKETS; i++) {
		file << i << "\t" << freq_buckets[i] << "\n";
	}
	file.close();
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
			for(uint32 k = 0; k < counts.size(); k++) {
				//if(counts[k] > 2) {
				//	repk += 1;
				//	break;
				//}
				window += std::to_string(counts[k]);
				window += " ";
			}
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

void kmer_stats(const char* refFname) {
	std::string fname(refFname);
	fname += std::string(".kmer_hist");
	std::ifstream file;
	file.open(fname.c_str(), std::ios::in | std::ios::binary);
	if (!file.is_open()) {
		printf("load_kmer_hist: Cannot open the hist file %s!\n", fname.c_str());
		exit(1);
	}

	std::string stats_fname(refFname);
	stats_fname += std::string("__kmer_stats.csv");
	std::ofstream stats_file;
	stats_file.open(stats_fname.c_str(), std::ios::out);
	if (!stats_file.is_open()) {
		printf("store_ref_idx: Cannot open the KMER STATS file %s!\n", stats_fname.c_str());
		exit(1);
	}
	uint32 kmer, count;
	int map_size;
	file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
	stats_file << "Count\n";
	while(map_size >= 0) {
		file.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
		file.read(reinterpret_cast<char*>(&count), sizeof(count));
		stats_file << count  << "\n";
		map_size--;
	}
	file.close();
}

#define BUCKET_SIZE_THR_DEBUG 50000
void store_ref_index_stats(const char* refFname, const ref_t& ref, const index_params_t* params) {
	std::string fname(refFname);
	fname += std::string(".idx.");
	fname += std::string("h");
	fname += std::to_string(params->h);
	fname += std::string("_T");
	fname += std::to_string(params->n_tables);
	fname += std::string("_b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("_w");
	fname += std::to_string(params->ref_window_size);
	fname += std::string("_p");
	fname += std::to_string(params->n_buckets_pow2);
	fname += std::string("_k");
	fname += std::to_string(params->k);
	fname += std::string("_H");
	fname += std::to_string(params->max_count);
	fname += std::string("__stats");

	std::ofstream file;
	file.open(fname.c_str(), std::ios::out | std::ios::app);
	if (!file.is_open()) {
		printf("store_ref_idx: Cannot open the STATS IDX file %s!\n", fname.c_str());
		exit(1);
	}

	for(uint32 i = 0; i < params->n_tables; i++) {
		const buckets_t& buckets = ref.mutable_index.per_table_buckets[i];
		for(uint32 j = 0; j < params->n_buckets; j++) {
			const VectorSeqPos& bucket = buckets.buckets_data_vectors[j];
			uint32 size = bucket.size();
			uint32 len_avg = 0;
			for(uint32 k = 0; k < size; k++) {
				len_avg += bucket[k].len;
#if(DEBUG)
				if(size > BUCKET_SIZE_THR_DEBUG) {
					printf("T %d b %d size %d pos %u \n", i, j, size, bucket[k].pos);
					for(seq_t x = 0; x < params->ref_window_size; x++) {
						printf("%c", iupacChar[(int) ref.seq[bucket[k].pos+x]]);
					}
					printf("\n");
				}
#endif
			}
			if(size > 0) {
				len_avg = len_avg/size;
			}
			// table id, bucket id, size, average contig len
			file << i << "," << j << "," << size << "," << len_avg << "\n";
		}
	}
	file.close();
}
