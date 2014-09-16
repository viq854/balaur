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
int is_inform_ref_window(char* seq, int len) {
	char c = seq[0];
	int count = 0; 
	int countN = 0; 
	for(int i = 1; i < len; i++) {
		if(seq[i] == c) {
			count++; //return 1;
		} 
		if(seq[i] == 4) {
			countN++;
		}
		if(count > len/2) {
			return 0;
		}
		if(countN > 5) {
			return 0;
		}
	}
	return 1;
}

// compute and store the frequency of each kmer in the given sequence
void compute_kmer_counts(const char* seq, const seq_t seq_len, const index_params_t* params,
		uint32_t* hist) {
	for(seq_t j = 0; j <= (seq_len - params->k); j++) {
		if(params->hist_size == KMER_HIST_SIZE16) {
			uint16_t kmer;
			if(pack_16(&seq[j], params->k, &kmer) < 0) {
				continue;
			}
			hist[kmer]++;
		} else {
			uint32_t kmer;
			if(pack_32(&seq[j], params->k, &kmer) < 0) {
				continue;
			}
			hist[kmer]++;
		}
	}
}

uint32_t get_kmer_count(const char* kmer_seq, int kmer_len, const uint32_t* hist,
		const index_params_t* params) {
	if(params->hist_size == KMER_HIST_SIZE16) {
		uint16_t kmer;
		if(pack_16(kmer_seq, kmer_len, &kmer) < 0) {
			return 0;
		}
		return hist[kmer];
	} else {
		uint32_t kmer;
		if(pack_32(kmer_seq, kmer_len, &kmer) < 0) {
			return 0;
		}
		return hist[kmer];
	}
}

// returns the weight of the given kmer
// 0 if the kmer should be ignored
uint32_t get_kmer_weight(const char* kmer_seq, int kmer_len,
		const uint32_t* read_hist, const uint32_t* ref_hist,
		const uint8_t is_ref,
		const index_params_t* params) {

	//const uint32_t ref_count = get_kmer_count(kmer_seq, kmer_len, ref_hist, params);
	//const uint32_t reads_count = get_kmer_count(kmer_seq, kmer_len, read_hist, params);

	// filter out kmers if:
	// 1. count is too low and kmer does not occur in the reference
	// 2. count is too high
	//if(/*(max_count == 0 && min_count < params->min_count) ||*/ (ref_count > params->max_count)) {
		//return 0;
	//}
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
hash_t simhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
		const uint32_t* reads_hist, const uint32_t* ref_hist,
		const index_params_t* params, const uint8_t is_ref) {

	int v[SIMHASH_BITLEN] = { 0 };

	// generate the kmers, hash them, and add the hash to V
	if(params->kmer_type == SPARSE) {
		char* kmer = (char*) malloc(params->k*sizeof(char));
		for(uint32_t i = 0; i < params->m; i++) {
			const uint32_t* ids = &params->sparse_kmers[i*params->k];
			for(uint32_t j = 0; j < params->k; j++) {
				kmer[j] = seq[seq_offset + ids[j]];
			}
			if(get_kmer_weight(kmer, params->k, reads_hist, ref_hist, is_ref, params) == 0) {
				continue;
			}
			hash_t kmer_hash = CityHash64(kmer, params->k);
			add_kmer_hash_bits(v, kmer_hash);
		}
	} else {
		for(uint32_t i = 0; i <= (seq_len - params->k); i += params->kmer_dist) {
			if(get_kmer_weight(&seq[seq_offset + i], params->k, reads_hist, ref_hist, is_ref, params) == 0) {
				continue;
			}
			hash_t kmer_hash = CityHash64(&seq[seq_offset + i], params->k);
			add_kmer_hash_bits(v, kmer_hash);
		}
	}
	return generate_simhash_fp(v);
}

/////////////////////////
// --- LSH: minhash ---

hash_t minhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
			const uint32_t* reads_hist, const uint32_t* ref_hist,
			const index_params_t* params, const uint8_t is_ref) {
	hash_t fingerprint = 0;
	char* kmer = (char*) malloc(params->k*sizeof(char));
	// find and store the min for each hash function
	for(uint32_t h = 0; h < params->h; h++) {
		hash_t min = LLONG_MAX;
		if(params->kmer_type == SPARSE) {
			// construct the kmers, hash them, and keep the min (only lowest bit)
			for(int i = 0; i < params->m; i++) {
				uint32_t* ids = &params->sparse_kmers[i*params->k];
				for(int j = 0; j < params->k; j++) {
					kmer[j] = seq[seq_offset + ids[j]];
				}
				// hash the k-mer and compare to current min
				hash_t kmer_hash = CityHash64(kmer, params->k);
				kmer_hash ^= params->rand_hash_pads[h]; // xor with the random pad
				if(kmer_hash < min) {
					min = kmer_hash;
				}
			}
		}
		else {
			for(int i = 0; i <= (seq_len - params->k); i += params->kmer_dist) {
				const char* kmer = &seq[seq_offset + i];
				int weight = get_kmer_weight(kmer, params->k, ref_hist, reads_hist, is_ref, params);
				if(weight == 0) continue;
				hash_t kmer_hash = CityHash64(kmer, params->k);
				kmer_hash ^= params->rand_hash_pads[h]; // xor with the random pad
				if(kmer_hash < min) {
					min = kmer_hash;
				}
			}
		}
		// keep only the lowest bit of the min
		fingerprint |= (min & 1ULL) << (1*h);
	}
	free(kmer);
	return fingerprint;
}

/////////////////////////
// --- LSH: sampling ---

hash_t sampling(const char* seq, const seq_t seq_offset, const seq_t i, const index_params_t* params) {
	hash_t fingerprint = 0;
	const uint32_t* idxs = &params->sparse_kmers[i*params->k];
	for(uint32_t j = 0; j < params->k; j++) {
		const char c = seq[seq_offset + idxs[j]];
		fingerprint |= (c & 1ULL) << j; // 1st ls bit
	}
	return fingerprint;
}

hash_t sampling_hash(const char* seq, const seq_t seq_offset, const seq_t i, const index_params_t* params) {
	const uint32_t* idxs = &params->sparse_kmers[i*params->k];
	char* kmer = (char*) malloc(params->k*sizeof(char));
	for(uint32_t j = 0; j < params->k; j++) {
		kmer[j] = seq[seq_offset + idxs[j]];
	}
	return CityHash64(kmer, params->k);
}


