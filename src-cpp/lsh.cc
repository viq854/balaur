#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "limits.h"
#include "lsh.h"
#include "hash.h"


// returns the hamming distance between two 64-bit fingerprints
int hamming_dist(hash_t h1, hash_t h2) {
	return __builtin_popcountll(h1 ^ h2);
}

// returns the weight of the given kmer
// 0 if the kmer should be ignored
uint32_t get_kmer_weight(const char* kmer_seq, uint32 kmer_len,
		const marisa::Trie& ref_high_freq_hist,
		const marisa::Trie& reads_low_freq_hist,
		const index_params_t* params) {

	for (uint32 k = 0; k < kmer_len; k++) {
		if(kmer_seq[k] == BASE_IGNORE) {
			return 0; // contains ambiguous bases
		}
	}

	// lookup high frequency kmer trie
	marisa::Agent agent;
	agent.set_query(kmer_seq, kmer_len);
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
// --- LSH: minhash ---

bool minhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
			const marisa::Trie& ref_freq_kmer_trie,
			const VectorBool& ref_freq_kmer_bitmask,
			const marisa::Trie& reads_hist,
			const index_params_t* params, CyclicHash* kmer_hasher, const bool is_ref,
			VectorMinHash& min_hashes) {

//	kmer_hasher->hashvalue = 0;
//	for(uint32 i = 0; i < params->k; i++) {
//		unsigned char c = seq[seq_offset + i];
//		kmer_hasher->eat(c);
//	}

//	bool any_valid_kmers = false;
//	for(uint32 i = 0; i <= (seq_len - params->k); i++) {
//		// check if the kmer should be discarded
//		if(!get_kmer_weight(&seq[seq_offset + i], params->k, ref_freq_kmer_trie, reads_hist, params)) {
//			//minhash_t kmer_hash = kmer_hasher->hashvalue;
//			minhash_t kmer_hash = CityHash32(&seq[seq_offset + i], params->k);
//			for(uint32_t h = 0; h < params->h; h++) { // update the min values
//				const rand_hash_function_t* f = &params->minhash_functions[h];
//				minhash_t min = f->apply(kmer_hash);
//				if(min < min_hashes[h] || i == 0) {
//					min_hashes[h] = min;
//				}
//			}
//			//any_valid_kmers = true;
//		}
//		// roll the hash
////		if(i < seq_len - params->k) {
////			unsigned char c_out = seq[seq_offset + i];
////			unsigned char c_in = seq[seq_offset + i + params->k];
////			kmer_hasher->update(c_out, c_in);
////		}
//	}
//	return any_valid_kmers;
	bool any_valid_kmers = false;
	for(uint32 i = 0; i <= (seq_len - params->k); i++) {
		// check if the kmer should be discarded
		char weight = 0;
		if(is_ref) {
			weight = !ref_freq_kmer_bitmask[seq_offset + i];
		} else {
			weight = get_kmer_weight(&seq[seq_offset + i], params->k, ref_freq_kmer_trie, reads_hist, params);
		}
		if(weight != 0) {
			//minhash_t kmer_hash = kmer_hasher->hashvalue;
			minhash_t kmer_hash = CityHash32(&seq[seq_offset + i], params->k);
			for(uint32_t h = 0; h < params->h; h++) { // update the min values
				const rand_hash_function_t* f = &params->minhash_functions[h];
				minhash_t min = f->apply(kmer_hash);
				if(min < min_hashes[h] || !any_valid_kmers) {
					min_hashes[h] = min;
				}
			}
		} else {
			continue;
		}

		// roll the hash
//		if(i < seq_len - params->k) {
//			unsigned char c_out = seq[seq_offset + i];
//			unsigned char c_in = seq[seq_offset + i + params->k];
//			kmer_hasher->update(c_out, c_in);
//		}
		any_valid_kmers = true;
	}

	if(!any_valid_kmers) {
		std::fill(min_hashes.begin(), min_hashes.end(), UINT_MAX);
	}
	return any_valid_kmers;
}

// avoid redundant computations
// reference-only
bool minhash_rolling_init(const char* seq, const seq_t ref_offset, const seq_t seq_len,
					minhash_matrix_t& rolling_minhash_matrix,
					const VectorBool& ref_freq_kmer_bitmask,
					const index_params_t* params, CyclicHash* kmer_hasher,
					VectorMinHash& min_hashes) {

	// initialize the rolling matrix
	rolling_minhash_matrix.h_minhash_cols.resize(seq_len - params->k + 1);
	for(uint32 pos = 0; pos < seq_len - params->k + 1; pos++) {
		rolling_minhash_matrix.h_minhash_cols[pos].resize(params->h);
	}
	rolling_minhash_matrix.oldest_col_index = 0;

	bool any_valid_kmers = false;
	for(uint32 i = 0; i < seq_len - params->k + 1; i++) {
		if(!ref_freq_kmer_bitmask[ref_offset + i]) { // check if the kmer should be discarded
			minhash_t kmer_hash = CityHash32(&seq[ref_offset + i], params->k);
			for(uint32_t h = 0; h < params->h; h++) { // update the min values
				const rand_hash_function_t* f = &params->minhash_functions[h];
				minhash_t min = f->apply(kmer_hash);
				rolling_minhash_matrix.h_minhash_cols[i][h] = min;
				if(min < min_hashes[h] || !any_valid_kmers) {
					min_hashes[h] = min;
				}
			}
			any_valid_kmers = true;
		} else {
			for(uint32_t h = 0; h < params->h; h++) {
				rolling_minhash_matrix.h_minhash_cols[i][h] = UINT_MAX;
			}
		}
	}

	if(!any_valid_kmers) {
		std::fill(min_hashes.begin(), min_hashes.end(), UINT_MAX);
	}

	return any_valid_kmers;
}

bool minhash_rolling(const char* seq, const seq_t ref_offset, const seq_t seq_len,
					minhash_matrix_t& rolling_minhash_matrix,
					const VectorBool& ref_freq_kmer_bitmask,
					const index_params_t* params, CyclicHash* kmer_hasher,
					VectorMinHash& min_hashes) {

	minhash_t new_kmer_hash = 0;
	bool new_kmer_hash_valid = false;
	seq_t last_kmer_pos = ref_offset + seq_len - params->k;
	if(!ref_freq_kmer_bitmask[last_kmer_pos]) { // check if the kmer should be discarded
		new_kmer_hash_valid = true;
		new_kmer_hash = CityHash32(&seq[last_kmer_pos], params->k);
	}
	bool any_valid_kmers = false;
	for(uint32 h = 0; h < params->h; h++) {
		minhash_t min_h = UINT_MAX;
		if(new_kmer_hash_valid) {
			const rand_hash_function_t* f = &params->minhash_functions[h];
			min_h = f->apply(new_kmer_hash);
		}
		// if less than current min, update min, populate oldest column
		if(min_h < min_hashes[h]) {
			min_hashes[h] = min_h;
			rolling_minhash_matrix.h_minhash_cols[rolling_minhash_matrix.oldest_col_index][h] = min_h;
		} else if(rolling_minhash_matrix.h_minhash_cols[rolling_minhash_matrix.oldest_col_index][h] != min_hashes[h]) {
			// if the minimum doesn't come from the old column, no need to recompute
			rolling_minhash_matrix.h_minhash_cols[rolling_minhash_matrix.oldest_col_index][h] = min_h;
		} else {
			// need to recompute the minimum
			min_hashes[h] = UINT_MAX;
			rolling_minhash_matrix.h_minhash_cols[rolling_minhash_matrix.oldest_col_index][h] = min_h;
			for(uint32 i = 0; i < seq_len - params->k + 1; i++) {
				if(rolling_minhash_matrix.h_minhash_cols[i][h] < min_hashes[h]) {
					min_hashes[h] = rolling_minhash_matrix.h_minhash_cols[i][h];
				}
			}
		}
		if(min_hashes[h] != UINT_MAX) {
			any_valid_kmers = true;
		}
	}
	rolling_minhash_matrix.oldest_col_index = (rolling_minhash_matrix.oldest_col_index + 1) % (seq_len - params->k + 1);
	return any_valid_kmers;
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


