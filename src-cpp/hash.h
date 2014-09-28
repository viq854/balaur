#ifndef HASH_H_
#define HASH_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include "types.h"
#include "io.h"
#include "index.h"
#include "city.h"

#define SIMHASH_BITLEN 64
#define MINHASH_BUCKET_SIZE 8

// hashing schemes
hash_t simhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
		const MapKmerCounts& ref_hist, const MapKmerCounts& reads_hist,
		const index_params_t* params, const uint8_t is_ref);
bool minhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
		const marisa::Trie& ref_freq_kmer_trie,
		const VectorBool& ref_freq_kmers_bitmask,
		const marisa::Trie& reads_hist,
		const index_params_t* params, CyclicHash* hasher, const bool is_ref,
		VectorMinHash& min_hashes);
bool minhash_rolling_init(ref_t& ref, const seq_t ref_offset, const seq_t seq_len,
		const VectorBool& ref_freq_kmer_bitmask,
		const index_params_t* params, CyclicHash* kmer_hasher,
		VectorMinHash& min_hashes);
bool minhash_rolling(ref_t& ref, const seq_t ref_offset, const seq_t seq_len,
		const VectorBool& ref_freq_kmer_bitmask,
		const index_params_t* params, CyclicHash* kmer_hasher,
		VectorMinHash& min_hashes);
hash_t sampling(const char* seq, const seq_t seq_offset, const seq_t i,
		const index_params_t* params);
hash_t sampling_hash(const char* seq, const seq_t seq_offset, const seq_t i,
		const index_params_t* params);

// kmer histogram
void compute_kmer_counts(const char* seq, const seq_t seq_len,
		const index_params_t* params, MapKmerCounts& hist);
void find_high_freq_kmers(const MapKmerCounts& hist, MapKmerCounts& high_freq_hist,
		const index_params_t* params);
void find_low_freq_kmers(const MapKmerCounts& hist, MapKmerCounts& low_freq_hist,
		const index_params_t* params);

// utils
int is_inform_ref_window(const char* seq, const uint32_t len);
int hamming_dist(hash_t h1, hash_t h2);

void cityhash(read_t* r);

void hashlittle2( 
  const void *key,       /* the key to hash */
  size_t      length,    /* length of the key */
  uint32_t   *pc,        /* IN: primary initval, OUT: primary hash */
  uint32_t   *pb); 

#endif /*HASH_H_*/
