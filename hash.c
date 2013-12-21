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

// compute the frequency of each kmer in the read set
void generate_reads_kmer_hist(reads_t* reads, index_params_t* params, int** hist) {
	// stores the counts of each kmer
	int* histogram = (int*) calloc(65536, sizeof(int));
	
	// add the contribution of each read
	for(int i = 0; i < reads->count; i++) {
		read_t* r = &reads->reads[i];
		for(int j = 0; j <= (r->len - params->k); j++) {
			// compress the kmer into 16 bits
			uint16_t kmer;
			if(pack_16(&r->seq[j], params->k, &kmer)) {
				continue;
			}
			//printf("compressed = %lx \n", kmer);
			// update the count
			histogram[kmer]++;
		}
	}
	*hist = histogram;
}

// checks if the window is informative or not
// e.g. non-informative windows: same character is repeated throughout the window (NN...N)
int is_valid_window(char* window, index_params_t* params) {
	char c = window[0];
	for(int i = 1; i < params->ref_window_size; i++) {
		if(window[i] != c) {
			return 1;
		}
	}
	return 0;
}

// compute the frequency of each kmer in the read set
void generate_ref_kmer_hist(ref_t* ref, index_params_t* params, int** hist) {
	// stores the counts of each kmer
	int* histogram = (int*) calloc(65536, sizeof(int));
	
	ref->num_windows = 0;
	seq_t max_num_windows = ref->len - params->ref_window_size + 1;
	ref->windows = (ref_win_t*) malloc(max_num_windows*sizeof(ref_win_t));
	seq_t i;
	for(i = 0; i < max_num_windows; i++) {
		if(is_valid_window(&ref->seq[i], params)) {
			ref_win_t* window = &ref->windows[ref->num_windows];
			window->pos = i;
			ref->num_windows++;	
		} else {
			// skip until at least 1 potential valid kmer
			i += params->k;
		}
	}
	ref->windows = (ref_win_t*) realloc(ref->windows, ref->num_windows*sizeof(ref_win_t));
	
	// compute the k-mers
	for(seq_t j = 0; j <= (ref->len - params->k); j++) {
		// compress the kmer into 16 bits
		uint16_t kmer;
		if(pack_16(&ref->seq[j], params->k, &kmer)) {
			continue;
		}
		// update the count
		histogram[kmer]++;
	}
	
	*hist = histogram;
}

// returns the weight of each read kmer
// 0 if the kmer should be ignored
int get_reads_kmer_weight(char* seq, int len, int* histogram, index_params_t* params) {
	uint16_t kmer;
	if(pack_16(seq, len, &kmer)) {
		return 0;
	}
	int count = histogram[kmer];
	//printf("count %d \n", histogram[kmer]);
	
	// filter out infrequent kmers
	if(count < params->min_count) {
		return 0;
	}
	return 1;
}

// returns the weight of each reference kmer
// 0 if the kmer should be ignored
int get_ref_kmer_weight(char* seq, int len, int* histogram, index_params_t* params) {
	uint16_t kmer;
	if(pack_16(seq, len, &kmer) < 0) {
		return 0;
	}
	int count = histogram[kmer];
	//printf("count %d \n", histogram[kmer]);
	
	// filter out infrequent kmers
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
void simhash_read(read_t* r, int* histogram, index_params_t* params) {	
	int v[SIMHASH_BITLEN] = { 0 };
	
	// find the read kmers, hash them, and add the hash to V
	int i;
	for(i = 0; i <= (r->len - params->k); i++) {
		int weight = get_reads_kmer_weight(&r->seq[i], params->k, histogram, params);
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
void simhash_ref(ref_t* ref, ref_win_t* window, int* histogram, index_params_t* params) {
	//printf("w = %llu \n", window->pos);
	int v[SIMHASH_BITLEN] = { 0 };
	// find the read kmers, hash them, and add the hash to V
	for(int i = 0; i <= (params->ref_window_size - params->k); i++) {
		int weight = get_ref_kmer_weight(&ref->seq[window->pos + i], params->k, histogram, params);
		if(weight == 0) continue;
		simhash_t kmer_hash = CityHash64(&ref->seq[window->pos + i], params->k);
		//printf("hash = %llx \n", kmer_hash);
		add_kmer_hash_bits(v, kmer_hash);
	}
	window->simhash = generate_simhash_fp(v);
}
