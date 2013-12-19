#ifndef HASH_H_
#define HASH_H_
#include "types.h"
#include "io.h"
#include "index.h"

// ---- simhash ----

#define SIMHASH_BITLEN 64
void simhash(read_t* r, int* histogram, index_params_t* params);
void generate_kmer_hist(reads_t* reads, index_params_t* params, int** histogram);
int hamming_dist(simhash_t h1, simhash_t h2);

void cityhash(read_t* r);

void hashlittle2( 
  const void *key,       /* the key to hash */
  size_t      length,    /* length of the key */
  uint32_t   *pc,        /* IN: primary initval, OUT: primary hash */
  uint32_t   *pb); 

#endif /*HASH_H_*/
