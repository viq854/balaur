#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "io.h"
#include "hash.h"
#include "city.h"


// for each bit position i in the kmer hash
// if hash[i] is 1: increment v[i]; otherwise, decrement v[i]
void add_kmer_hash_bits(int* v, simhash_t hash) {
	for(int b = 0; b < SIMHASH_BITLEN; b++) {
		if(((hash >> b) & 1) == 1) {
			v[b]++;
		} else {
			v[b]--;
		}
	}
}

// returns the simhash fingerprint
void generate_simhash_fp(read_t* r, int* v) {
	simhash_t simhash = 0;
	for (int b = 0; b < SIMHASH_BITLEN; b++) {
		if(v[b] >= 0) {
			simhash |= 1 << b;
			r->simhash_popc++;
		}
	}
	r->simhash = simhash;
	
}

// computes the Charikar simhash fingerprint
void simhash(read_t* r, index_params_t* params) {	
	int v[SIMHASH_BITLEN] = { 0 };
	
	// find the read kmers, hash them, and add the hash to V
	for(int i = 0; i <= (r->len - params->k); i++) {
		simhash_t kmer_hash = CityHash64(&r->seq[i], params->k);
		//printf("hash = %llx \n", kmer_hash);
		add_kmer_hash_bits(v, kmer_hash);
	}
	generate_simhash_fp(r, v);
}