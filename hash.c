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

// --- hamming distance utils ---

int hamming_dist(simhash_t h1, simhash_t h2) {
	return __builtin_popcountl(h1 ^ h2);
}

// --- k-mer weights ---

// checks if the sequence is informative or not
// e.g. non-informative seq: same character is repeated throughout the seq (NN...N)
int is_inform(char* seq, int len) {
	char c = seq[0];
	for(int i = 1; i < len; i++) {
		if(seq[i] != c) {
			return 1;
		}
	}
	return 0;
}

// compute the frequency of each kmer in the read set
void generate_reads_kmer_hist(reads_t* reads, index_params_t* params) {
	// stores the counts of each kmer
	reads->hist = (int*) calloc(params->hist_size, sizeof(int));
	
	// add the contribution of each read
	for(int i = 0; i < reads->count; i++) {
		read_t* r = &reads->reads[i];
		for(int j = 0; j <= (r->len - params->k); j++) {
			if(params->hist_size == KMER_HIST_SIZE16) {
				// compress the kmer into 16 bits
				uint16_t kmer;
				if(pack_16(&r->seq[j], params->k, &kmer) < 0) {
					continue;
				}
				// update the count
				reads->hist[kmer]++;
			} else {
				// compress the kmer into 16 bits
				uint32_t kmer;
				if(pack_32(&r->seq[j], params->k, &kmer) < 0) {
					continue;
				}
				// update the count
				reads->hist[kmer]++;
			}
		}
	}
}

// compute the frequency of each kmer in the reference
void generate_ref_kmer_hist(ref_t* ref, index_params_t* params) {
	// stores the counts of each kmer
	ref->hist = (int*) calloc(params->hist_size, sizeof(int));
	ref->hist_size = params->hist_size;
	// compute the k-mers
	for(seq_t j = 0; j <= (ref->len - params->k); j++) {
		if(params->hist_size == KMER_HIST_SIZE16) {
			uint16_t kmer;
			int r = pack_16(&ref->seq[j], params->k, &kmer);
			// valid window
			if((j < ref->len - params->ref_window_size) && (is_inform(&ref->seq[j], params->ref_window_size))) {
				if(r >= 0) {
					ref->hist[kmer]++;
				}
			} else { // repetitive window
				if(r >= 0) {
					//not a stretch of N's
					ref->hist[kmer] += params->ref_window_size - params->k;	
				}
				j += params->ref_window_size - params->k;
			}
		} else {
			uint32_t kmer;
			int r = pack_32(&ref->seq[j], params->k, &kmer);
			// valid window
			if((j < ref->len - params->ref_window_size) && (is_inform(&ref->seq[j], params->ref_window_size))) {
				if(r >= 0) {
					ref->hist[kmer]++;
				}
			} else { // repetitive window
				if(r >= 0) {
					//not a stretch of N's
					ref->hist[kmer] += params->ref_window_size - params->k;	
				}
				j += params->ref_window_size - params->k;
			}
		}
	}
}

// returns the weight of each read kmer
// 0 if the kmer should be ignored
int get_reads_kmer_weight(char* seq, int len, int* reads_hist, int* ref_hist, index_params_t* params) {
	int min_count, max_count;
	if(params->hist_size == KMER_HIST_SIZE16) {
		uint16_t kmer;
		if(pack_16(seq, len, &kmer) < 0) {
			return 0;
		}
		min_count = reads_hist[kmer];
		max_count = ref_hist[kmer];
	} else {
		uint32_t kmer;
		if(pack_32(seq, len, &kmer) < 0) {
			return 0;
		}
		min_count = reads_hist[kmer];
		max_count = ref_hist[kmer];
	}
	
	//printf("count %d %llu \n", reads_hist[kmer], params->min_count);
	//printf("count %d %llu \n", ref_hist[kmer], params->max_count);
	
	// filter out kmers if:
	// 1. count is too low and kmer does not occur in the reference
	// 2. count is too high
	if((max_count == 0 && min_count < params->min_count) || (max_count > params->max_count)) {
		return 0;
	}
	return 1;
}

// returns the weight of each reference kmer
// 0 if the kmer should be ignored
int get_ref_kmer_weight(char* seq, int len, int* hist, index_params_t* params) {
	int count;
	if(params->hist_size == KMER_HIST_SIZE16) {
		uint16_t kmer;
		if(pack_16(seq, len, &kmer) < 0) {
			return 0;
		}
		count = hist[kmer];
	} else {
		uint32_t kmer;
		if(pack_32(seq, len, &kmer) < 0) {
			return 0;
		}
		count = hist[kmer];
	}
	//printf("count %d \n", hist[kmer]);
	
	// filter out kmers that are too frequent
	if(count > params->max_count) {
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
simhash_t generate_simhash_fp(int* v) {
	simhash_t simhash = 0;
	for (int b = 0; b < SIMHASH_BITLEN; b++) {
		if(v[b] >= 0) {
			simhash |= (1ULL << b);
			//r->simhash_popc++;
		}
	}
	return simhash;
}

// computes the Charikar simhash fingerprint
void simhash_read(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params) {	
	int v[SIMHASH_BITLEN] = { 0 };
	
	// find the read kmers, hash them, and add the hash to V
	int i;
	for(i = 0; i <= (r->len - params->k); i++) {
		int weight = get_reads_kmer_weight(&r->seq[i], params->k, reads_hist, ref_hist, params);
		if(weight == 0) continue;
		simhash_t kmer_hash = CityHash64(&r->seq[i], params->k);
		//printf("hash = %llx \n", kmer_hash);
		add_kmer_hash_bits(v, kmer_hash);
	}
	//add_kmer_hash_bits(v, CityHash64(&r->seq[i], (params->k - 1)));
	r->simhash = generate_simhash_fp(v);
}

// --- regular hashing ---
void cityhash(read_t* r) {
	r->simhash = CityHash64(r->seq, r->len);
}

// compute the simhash of a reference window
void simhash_ref(ref_t* ref, ref_win_t* window, index_params_t* params) {
	//if(window->pos == 26738515) {
		 //for(seq_t i = 0; i < 100; i++) {
			 //printf("%c", iupacChar[(int) ref->seq[window->pos + i]]);
		 //} 
		 //printf("\n");
	//}
	
	//printf("w = %llu \n", window->pos);
	int v[SIMHASH_BITLEN] = { 0 };
	// find the read kmers, hash them, and add the hash to V
	for(int i = 0; i <= (params->ref_window_size - params->k); i++) {
		int weight = get_ref_kmer_weight(&ref->seq[window->pos + i], params->k, ref->hist, params);
		if(weight == 0) continue;
		simhash_t kmer_hash = CityHash64(&ref->seq[window->pos + i], params->k);
		//printf("hash = %llx \n", kmer_hash);
		add_kmer_hash_bits(v, kmer_hash);
	}
	window->simhash = generate_simhash_fp(v);
	
	//if(window->pos == 26738515) {
		//printf("window hash %llx \n", window->simhash); 
	//}
}
