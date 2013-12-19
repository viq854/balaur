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

// --- k-mer weights ---

// compute the frequency of each kmer
void generate_kmer_hist(reads_t* reads, index_params_t* params, int** hist) {
	// stores the counts of each kmer
	int* histogram = (int*) calloc(65536, sizeof(int));
	
	// add the contribution of each read
	for(int i = 0; i < reads->count; i++) {
		read_t* r = &reads->reads[i];
		for(int j = 0; j <= (r->len - params->k); j++) {
			// compress the kmer into 16 bits
			uint16_t kmer;
			int e = pack_16(&r->seq[j], params->k, &kmer);
			if(e < 0) {
				continue; // ignore kmers that contain Ns
			}
			//printf("compressed = %lx \n", kmer);
			
			// update the count
			histogram[kmer]++;
		}
	}
	*hist = histogram;
}

// returns the weight of each kmer
// 0 if the kmer should be ignored
int get_kmer_weight(char* seq, int len, int* histogram, index_params_t* params) {
	uint16_t kmer;
	int e = pack_16(seq, len, &kmer);
	if(e < 0) {
		return 0; // ignore kmers that contain Ns
	}
	int count = histogram[kmer];
	//printf("count %d \n", histogram[kmer]);
	
	// filter out infrequent kmers
	if(count < params->min_count) {
		return 0;
	}
	return 1;
}

// --- LSH: simhash ---


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
			simhash |= (1ULL << b);
			r->simhash_popc++;
		}
	}
	r->simhash = simhash;
	
}

// computes the Charikar simhash fingerprint
void simhash(read_t* r, int* histogram, index_params_t* params) {	
	int v[SIMHASH_BITLEN] = { 0 };
	
	// find the read kmers, hash them, and add the hash to V
	int i;
	for(i = 0; i <= (r->len - params->k); i++) {
		int weight = get_kmer_weight(&r->seq[i], params->k, histogram, params);
		if(weight == 0) continue;
		simhash_t kmer_hash = CityHash64(&r->seq[i], params->k);
		//printf("hash = %llx \n", kmer_hash);
		add_kmer_hash_bits(v, kmer_hash);
	}
	//add_kmer_hash_bits(v, CityHash64(&r->seq[i], (params->k - 1)));
	generate_simhash_fp(r, v);
}

// --- regular hashing ---
void cityhash(read_t* r) {
	r->simhash = CityHash64(r->seq, r->len);
}