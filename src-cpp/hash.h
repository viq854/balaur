#ifndef HASH_H_
#define HASH_H_

#include <vector>
#include "types.h"
#include "io.h"
#include "index.h"

#define SIMHASH_BITLEN 64
#define MINHASH_BUCKET_SIZE 8

// hashing schemes
hash_t simhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
		const MapKmerCounts& reads_hist, const MapKmerCounts& ref_hist,
		const index_params_t* params, const uint8_t is_ref);
hash_t minhash(const char* seq, const seq_t seq_offset, const seq_t seq_len,
		const MapKmerCounts& reads_hist, const MapKmerCounts& ref_hist,
		const index_params_t* params, const uint8_t is_ref, VectorMinHash& min_hashes);
hash_t sampling(const char* seq, const seq_t seq_offset, const seq_t i,
		const index_params_t* params);
hash_t sampling_hash(const char* seq, const seq_t seq_offset, const seq_t i,
		const index_params_t* params);

// kmer histogram
void compute_kmer_counts(const char* seq, const seq_t seq_len,
		const index_params_t* params, MapKmerCounts& hist);

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
