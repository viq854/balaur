#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "limits.h"
#include "io.h"
#include "hash.h"
#include "city.h"

// --- hamming distance utils ---

int hamming_dist(hash_t h1, hash_t h2) {
	return __builtin_popcountll(h1 ^ h2);
}

// --- k-mer weights ---

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

uint32_t get_kmer_count(const char* kmer_seq, int kmer_len, const MapKmerCounts& hist,
		const index_params_t* params) {
	uint32_t kmer;
	if(pack_32(kmer_seq, kmer_len, &kmer) < 0) {
		return 0;
	}
	MapKmerCounts::const_iterator v;
	if((v = hist.find(kmer)) != hist.end()) {
		return v->second;
	} else {
		return 0;
	}
}

void find_high_freq_kmers(const MapKmerCounts& hist, MapKmerCounts& high_freq_hist,
		const index_params_t* params) {

	seq_t total_counts = 0;
	seq_t max_count = 0;
	seq_t filt_count = 0;
	for(MapKmerCounts::const_iterator it = hist.begin(); it != hist.end(); ++it) {
		seq_t count = it->second;
		if(count > params->max_count) {
			high_freq_hist[it->first]++;
			filt_count++;
		}

		total_counts += count;
		if(count > max_count) {
			max_count = count;
		}
	}

	float avg_count = (float) total_counts/hist.size();
	printf("Total count %u %zu \n", total_counts, hist.size());
	printf("Avg count %.4f \n", avg_count);
	printf("Max count %u \n", max_count);
	printf("Count above cutoff %u \n", filt_count);
}

void find_low_freq_kmers(const MapKmerCounts& hist, MapKmerCounts& low_freq_hist,
		const index_params_t* params) {

	for(MapKmerCounts::const_iterator it = hist.begin(); it != hist.end(); ++it) {
		seq_t count = it->second;
		if(count < params->min_count) {
			low_freq_hist[it->first]++;
		}
	}
}

// returns true if the kmer occurs with high frequency in the reference
bool contains_kmer(const uint32_t kmer, const MapKmerCounts& freq_hist) {

	MapKmerCounts::const_iterator v;
	if((v = freq_hist.find(kmer)) != freq_hist.end()) {
		return true;
	}
	return false;
}

// returns the weight of the given kmer
// 0 if the kmer should be ignored
uint32_t get_kmer_weight(const char* kmer_seq, uint32 kmer_len,
		const marisa::Trie& ref_high_freq_hist,
		const marisa::Trie& reads_low_freq_hist,
		const uint8_t is_ref,
		const index_params_t* params) {

	for (uint32 k = 0; k < kmer_len; k++) {
		if(kmer_seq[k] == BASE_IGNORE) {
			return 0; // contains ambiguous bases
		}
	}

	// lookup high frequency kmer trie
	marisa::Agent agent;
	agent.set_query(kmer_seq);
	if(ref_high_freq_hist.lookup(agent)) {
		return 0;
	}

//	if(!is_ref) {
//		if(contains_kmer(kmer, reads_low_freq_hist)) {
//			return 0;
//		}
//	}
	return 1;
}

/////////////////////////
// --- LSH: simhash ---


// for each bit position i in the kmer hash
// if hash[i] is 1: increment v[i]; otherwise, decrement v[i]
void add_kmer_hash_bits(int* v, hash_t hash) {
	for(int b = 0; b < SIMHASH_BITLEN; b++) {
		if(((hash >> b) & 1) == 1) {
			v[b]++;
		} else {
			v[b]--;
		}
	}
}

// computes the simhash fingerprint
hash_t generate_simhash_fp(int* v) {
	hash_t simhash = 0;
	for (int b = 0; b < SIMHASH_BITLEN; b++) {
		if(v[b] >= 0) {
			simhash |= (1ULL << b);
		}
	}
	return simhash;
}

// computes the simhash fingerprint of the given sequence
// using the specified kmer generation scheme
//hash_t simhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
//		const MapKmerCounts& ref_hist, const MapKmerCounts& reads_hist,
//		const index_params_t* params, const uint8_t is_ref) {
//
//	int v[SIMHASH_BITLEN] = { 0 };
//
//	// generate the kmers, hash them, and add the hash to V
//	if(params->kmer_type == SPARSE) {
//		char* kmer = (char*) malloc(params->k*sizeof(char));
//		for(uint32_t i = 0; i < params->m; i++) {
//			const uint32_t* ids = &params->sparse_kmers[i*params->k];
//			for(uint32_t j = 0; j < params->k; j++) {
//				kmer[j] = seq[seq_offset + ids[j]];
//			}
//			if(get_kmer_weight(kmer, params->k, ref_hist, reads_hist, is_ref, params) == 0) {
//				continue;
//			}
//			hash_t kmer_hash = CityHash64(kmer, params->k);
//			add_kmer_hash_bits(v, kmer_hash);
//		}
//	} else {
//		for(uint32_t i = 0; i <= (seq_len - params->k); i += params->kmer_dist) {
//			if(get_kmer_weight(&seq[seq_offset + i], params->k, ref_hist, reads_hist, is_ref, params) == 0) {
//				continue;
//			}
//			hash_t kmer_hash = CityHash64(&seq[seq_offset + i], params->k);
//			add_kmer_hash_bits(v, kmer_hash);
//		}
//	}
//	return generate_simhash_fp(v);
//}

/////////////////////////
// --- LSH: minhash ---

bool minhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
			const marisa::Trie& ref_hist, const marisa::Trie& reads_hist,
			const index_params_t* params, CyclicHash* kmer_hasher, const uint8_t is_ref,
			VectorMinHash& min_hashes) {

	kmer_hasher->hashvalue = 0;
	for(uint32 i = 0; i < params->k; i++) {
		unsigned char c = seq[seq_offset + i];
		kmer_hasher->eat(c);
	}

	std::fill(min_hashes.begin(), min_hashes.end(), UINT32_MAX);
	for(uint32 i = 0; i <= (seq_len - params->k); i++) {
		// check if the kmer should be discarded
		if(get_kmer_weight(&seq[seq_offset + i], params->k, ref_hist, reads_hist, is_ref, params)) continue;
		minhash_t kmer_hash = kmer_hasher->hashvalue;

		// update the mins
		for(uint32_t h = 0; h < params->h; h++) {
			const rand_hash_function_t* f = &params->minhash_functions[h];
			minhash_t min = f->apply(kmer_hash);
			if(min < min_hashes[h]) {
				min_hashes[h] = min;
			}
		}

		// roll the hash
		if(i < seq_len - params->k) {
			unsigned char c_out = seq[seq_offset + i];
			unsigned char c_in = seq[seq_offset + i + params->k];
			kmer_hasher->update(c_out, c_in);
		}
	}
	return !(min_hashes[0] == UINT32_MAX);
}


/////////////////////////
// --- LSH: sampling ---

//hash_t sampling(const char* seq, const seq_t seq_offset, const seq_t i, const index_params_t* params) {
//	hash_t fingerprint = 0;
//	const uint32_t* idxs = &params->sparse_kmers[i*params->k];
//	for(uint32_t j = 0; j < params->k; j++) {
//		const char c = seq[seq_offset + idxs[j]];
//		fingerprint |= (c & 1ULL) << j; // 1st ls bit
//	}
//	return fingerprint;
//}
//
//hash_t sampling_hash(const char* seq, const seq_t seq_offset, const seq_t i, const index_params_t* params) {
//	const uint32_t* idxs = &params->sparse_kmers[i*params->k];
//	char* kmer = (char*) malloc(params->k*sizeof(char));
//	for(uint32_t j = 0; j < params->k; j++) {
//		kmer[j] = seq[seq_offset + idxs[j]];
//	}
//	return CityHash64(kmer, params->k);
//}


