#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "limits.h"
#include "io.h"
#include "hash.h"
#include "city.h"

// --- hamming distance utils ---

int hamming_dist(simhash_t h1, simhash_t h2) {
	return __builtin_popcountll(h1 ^ h2);
}

// --- k-mer weights ---

// checks if the sequence is informative or not
// e.g. non-informative seq: same character is repeated throughout the seq (NN...N)
int is_inform(char* seq, int len) {
	char c = seq[0];
	int count = 0; 
	int countN = 0; 
	for(int i = 1; i < len; i++) {
		if(seq[i] == c) {
			count++; //return 1;
		} 
		if(seq[i] == 4) {
			countN++;
		}
		if(count > len/2) {
			return 0;
		}
		if(countN > 5) {
			return 0;
		}
	}
	return 1;
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
				// compress the kmer into 32 bits
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
			if(pack_16(&ref->seq[j], params->k, &kmer)) {
				ref->hist[kmer]++;
			}
		} else {
			uint32_t kmer;
			if(pack_32(&ref->seq[j], params->k, &kmer)) {
				ref->hist[kmer]++;
			}
		}
	}
}

// compute the frequency of each kmer in the reference
void generate_ref_kmer_hist_sparse(ref_t* ref, index_params_t* params) {
	// stores the counts of each kmer
	ref->hist = (int*) calloc(params->hist_size, sizeof(int));
	ref->hist_size = params->hist_size;
	// compute the k-mers
	char* kmer_seq = (char*) malloc(params->k*sizeof(char));
	
	for(int pos = 0; pos < ref->len; pos += params->ref_window_size) {
		for(int i = 0; i < params->m; i++) {
			int* ids = &params->sparse_kmers[i*params->k]; 
			for(int j = 0; j < params->k; j++) {
				kmer_seq[j] = ref->seq[pos + ids[j]];
			}
			if(params->hist_size == KMER_HIST_SIZE16) {
				uint16_t kmer;
				if(pack_16(kmer_seq, params->k, &kmer) >= 0) {
					ref->hist[kmer]++;
				}
			} else {
				uint32_t kmer;
				if(pack_32(kmer_seq, params->k, &kmer) >= 0) {
					ref->hist[kmer]++;
				}
			}
		}
	}
	free(kmer_seq);
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
	if(/*(max_count == 0 && min_count < params->min_count) ||*/ (max_count > params->max_count)) {
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

// computes the simhash fingerprint using overlapping k-mers
void simhash_read_ovp(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params) {	
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

// computes the simhash fingerprint using non-overlapping k-mers
void simhash_read_novp(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params) {	
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

// compute the simhash fingerprint using sparse k-mers
void simhash_read_sparse(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params) {
	int v[SIMHASH_BITLEN] = { 0 };
	// find the kmers, hash them, and add the hash to V
	char* kmer = (char*) malloc(params->k*sizeof(char));
	for(int i = 0; i <params->m; i++) {
		int* ids = &params->sparse_kmers[i*params->k]; 
		for(int j = 0; j < params->k; j++) {
			kmer[j] = r->seq[ids[j]];
		}
		int weight = get_reads_kmer_weight(kmer, params->k, reads_hist, ref_hist, params);
		if(weight == 0) continue;
		simhash_t kmer_hash = CityHash64(kmer, params->k);
		add_kmer_hash_bits(v, kmer_hash);
	}
	r->simhash = generate_simhash_fp(v);
}

// compute the simhash of a reference window using overlapping k-mers
void simhash_ref_ovp(ref_t* ref, ref_win_t* window, index_params_t* params) {
	int v[SIMHASH_BITLEN] = { 0 };
	// find the kmers, hash them, and add the hash to V
	for(int i = 0; i <= (params->ref_window_size - params->k); i++) {
		int weight = get_ref_kmer_weight(&ref->seq[window->pos + i], params->k, ref->hist, params);
		if(weight == 0) continue;
		simhash_t kmer_hash = CityHash64(&ref->seq[window->pos + i], params->k);
		//printf("hash = %llx \n", kmer_hash);
		add_kmer_hash_bits(v, kmer_hash);
	}
	window->simhash = generate_simhash_fp(v);
}

// compute the simhash of a reference window using non-overlapping k-mers
void simhash_ref_novp(ref_t* ref, ref_win_t* window, index_params_t* params) {
	int v[SIMHASH_BITLEN] = { 0 };
	// find the kmers, hash them, and add the hash to V
	for(int i = 0; i <= (params->ref_window_size - params->k); i += params->k) {
		int weight = get_ref_kmer_weight(&ref->seq[window->pos + i], params->k, ref->hist, params);
		if(weight == 0) continue;
		simhash_t kmer_hash = CityHash64(&ref->seq[window->pos + i], params->k);
		add_kmer_hash_bits(v, kmer_hash);
	}
	window->simhash = generate_simhash_fp(v);
}

// compute the simhash of a reference window using sparse k-mers
void simhash_ref_sparse(ref_t* ref, ref_win_t* window, index_params_t* params) {
	int v[SIMHASH_BITLEN] = { 0 };
	// find the kmers, hash them, and add the hash to V
	char* kmer = (char*) malloc(params->k*sizeof(char));
	for(int i = 0; i <params->m; i++) {
		int* ids = &params->sparse_kmers[i*params->k]; 
		for(int j = 0; j < params->k; j++) {
			kmer[j] = ref->seq[window->pos + ids[j]];
			//printf("%d", kmer[j]);
		}
		//printf("\n");
		int weight = get_ref_kmer_weight(kmer, params->k, ref->hist, params);
		if(weight == 0) continue;
		simhash_t kmer_hash = CityHash64(kmer, params->k);
		//printf("kmer hash = %llx \n", kmer_hash);
		add_kmer_hash_bits(v, kmer_hash);
	}
	window->simhash = generate_simhash_fp(v);
	//printf("hash = %llx \n", window->simhash);
}

// --- LSH Sampling ---

void sampling_ref(ref_t* ref, ref_win_t* window, index_params_t* params, int i) {
	window->simhash = 0;
	int* idxs = &params->sparse_kmers[i*params->k]; 
	for(int j = 0; j < params->k; j++) {
		char c = ref->seq[window->pos + idxs[j]];
		window->simhash |= (c & 1ULL) << j; // 1st ls bit
	}
}

void sampling_read(read_t* r, index_params_t* params, int i) {
	r->simhash = 0;
	int* idxs = &params->sparse_kmers[i*params->k]; 
	for(int j = 0; j < params->k; j++) {
		char c = r->seq[idxs[j]];
		r->simhash |= (c & 1ULL) << j; // 1st ls bit
	}
}

void sampling_hash_ref(ref_t* ref, ref_win_t* window, index_params_t* params, int i) {
	int* idxs = &params->sparse_kmers[i*params->k]; 
	char* kmer = (char*) malloc(params->k*sizeof(char));
	for(int j = 0; j < params->k; j++) {
		kmer[j] = ref->seq[window->pos + idxs[j]];
	}
	window->simhash = CityHash64(kmer, params->k);
}

void sampling_hash_read(read_t* r, index_params_t* params, int i) {
	int* idxs = &params->sparse_kmers[i*params->k]; 
	char* kmer = (char*) malloc(params->k*sizeof(char));
	for(int j = 0; j < params->k; j++) {
		kmer[j] = r->seq[idxs[j]];
	}
	r->simhash = CityHash64(kmer, params->k);
}

// --- MIN Hash ---

void minhash_ref(ref_t* ref, ref_win_t* window, index_params_t* params) {
	window->simhash = 0;
	// find the kmers, hash them, and keep the min (only lowest bit)
	for(int j = 0; j < params->h; j++) {
		simhash_t min = LLONG_MAX; //INT_MAX; 
		// generate all non-overlapping windows
		for(int i = 0; i <= (params->ref_window_size - params->k); i += params->k) {
			// TODO generate sparse k-mers from each window
			int weight = get_ref_kmer_weight(&ref->seq[window->pos + i], params->k, ref->hist, params);
			if(weight == 0) continue;
			// hash the k-mer and compare to current min
			simhash_t kmer_hash = CityHash64(&ref->seq[window->pos + i], params->k);
			if(i > 0) {
				kmer_hash ^= params->rand_hash_pads[j]; // xor with the random pad 
			}
			if(kmer_hash < min) {
				min = kmer_hash;
			}
//			printf("kmer %llx \n", kmer_hash);
		}
		// keep only the lowest bit of the min
//		printf("max %llx min %llx bit %d \n", INT_MAX, min, (min & 1ULL));
		window->simhash |= (min & 1ULL) << (1*j); 
	}
//	printf("%llx \n", window->simhash);
}

void minhash_read(read_t* r, int* reads_hist, int* ref_hist, index_params_t* params) {
	r->simhash = 0;
	// find the kmers, hash them, and keep the min (only lowest bit)
	for(int j = 0; j < params->h; j++) {
		simhash_t min = LLONG_MAX; //INT_MAX; 
		// generate all non-overlapping windows
		for(int i = 0; i <= (r->len - params->k); i += params->k) {
			// TODO generate sparse k-mers from each window
			int weight = get_reads_kmer_weight(&r->seq[i], params->k, reads_hist, ref_hist, params);
			if(weight == 0) continue;
			// hash the k-mer and compare to current min
			simhash_t kmer_hash = CityHash64(&r->seq[i], params->k);
			if(i > 0) {
				kmer_hash ^= params->rand_hash_pads[j]; // xor with the random pad 
			}
			if(kmer_hash < min) {
				min = kmer_hash;
			}
		}
		// keep only the lowest bit of the min
		r->simhash |= (min & 1ULL) << (1*j); 
	}
}

