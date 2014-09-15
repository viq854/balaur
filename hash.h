#ifndef HASH_H_
#define HASH_H_
#include "types.h"
#include "io.h"
#include "index.h"

#define SIMHASH_BITLEN 64
#define MINHASH_BUCKET_SIZE 8

// simhash
void simhash_read_ovp(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params);
void simhash_read_novp(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params);
void simhash_read_sparse(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params);
void simhash_ref_ovp(ref_t* ref, ref_win_t* w, index_params_t* params);
void simhash_ref_novp(ref_t* ref, ref_win_t* w, index_params_t* params);
void simhash_ref_sparse(ref_t* ref, ref_win_t* w, index_params_t* params);

// sampling
void sampling_hash_read(read_t* r, index_params_t* params, int i);
void sampling_hash_ref(ref_t* ref, ref_win_t* window, index_params_t* params, int i);
void sampling_read(read_t* r, index_params_t* params, int i);
void sampling_ref(ref_t* ref, ref_win_t* window, index_params_t* params, int i);

// minhash
void minhash_ref(ref_t* ref, ref_win_t* window, index_params_t* params);
void minhash_read(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params);

// histograms
void generate_reads_kmer_hist(reads_t* reads, index_params_t* params);
void generate_reads_kmer_hist_sparse(reads_t* reads, index_params_t* params);
void generate_ref_kmer_hist(ref_t* ref, index_params_t* params);
void generate_ref_kmer_hist_sparse(ref_t* ref, index_params_t* params) ;

int is_inform(char* seq, int len);
int hamming_dist(hash_t h1, hash_t h2);

void cityhash(read_t* r);

void hashlittle2( 
  const void *key,       /* the key to hash */
  size_t      length,    /* length of the key */
  uint32_t   *pc,        /* IN: primary initval, OUT: primary hash */
  uint32_t   *pb); 

#endif /*HASH_H_*/
