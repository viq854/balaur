#ifndef LSH_H_
#define LSH_H_

#pragma once

#include "types.h"
#include "io.h"
#include "index.h"

#define SIMHASH_BITLEN 64

// LSH schemes

// rolling min-hash structure
struct minhash_matrix_t {
	std::vector<VectorMinHash> h_minhash_cols;
	uint32 oldest_col_index;
};

bool minhash(const char* seq, const seq_t seq_len,
		const VectorBool& ref_freq_kmer_bitmask,
		const MarisaTrie& ref_freq_kmer_trie,
		const MarisaTrie& reads_hist,
		const index_params_t* params, CyclicHash* hasher,
		VectorMinHash& min_hashes);
bool minhash_rolling_init(const char* seq, const seq_t ref_offset, const seq_t seq_len,
		minhash_matrix_t& rolling_minhash_matrix,
		const VectorBool& ref_freq_kmer_bitmask,
		const index_params_t* params, CyclicHash* kmer_hasher,
		VectorMinHash& min_hashes);
bool minhash_rolling(const char* seq, const seq_t ref_offset, const seq_t seq_len,
		minhash_matrix_t& rolling_minhash_matrix,
		const VectorBool& ref_freq_kmer_bitmask,
		const index_params_t* params, CyclicHash* kmer_hasher,
		VectorMinHash& min_hashes);

hash_t simhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
		const MapKmerCounts& ref_hist, const MapKmerCounts& reads_hist,
		const index_params_t* params, const uint8_t is_ref);
hash_t sampling(const char* seq, const seq_t seq_offset, const seq_t i,
		const index_params_t* params);
hash_t sampling_hash(const char* seq, const seq_t seq_offset, const seq_t i,
		const index_params_t* params);


#endif /* LSH_H_ */
