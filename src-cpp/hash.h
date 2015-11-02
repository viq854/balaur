#ifndef HASH_H_
#define HASH_H_

#pragma once
#include "types.h"
#include "city.h"
#include "mt64.h"
//#include "../rollinghashcpp/cyclichash.h"
struct CyclicHash{ CyclicHash(int x, int y){}; };

void sha1_hash(const uint8_t *message, uint32_t len, uint32_t hash[5]);

// universal hash function:
// single value: a*x mod n_buckets
// vector value: sum(a[i]*x[i]) mod n_buckets
struct rand_hash_function_t {
	const static uint32 w = 32; // bits per word
	minhash_t a; // odd a < 2^w
	std::vector<uint64> a_vec;
	minhash_t M; // n_butckets = 2^M

	// single int hashing
	rand_hash_function_t() {
		a = rand() + 1;
		M = w;
	}

	// vector hashing
	rand_hash_function_t(minhash_t _M, const uint32 vec_len) {
		a = 0;
		for(uint32 i = 0; i < vec_len; i++) {
			a_vec.push_back(genrand64_int64() + 1);
		}
		M = _M;
	}

	minhash_t apply(minhash_t x) const {
		return (minhash_t) a*x >> (w - M);
	}

	minhash_t apply_vector(const VectorMinHash& x, const VectorU32& indices, const uint32 vec_offset) const {
		uint64 s = 0;
		for(uint32 i = 0; i < a_vec.size(); i++) {
			s += a_vec[i]*x[indices[vec_offset + i]];
		}
		return (minhash_t) s;
		//return (minhash_t) s >> (w - M);
	}

	minhash_t bucket_hash(const minhash_t vector_prod) const {
		return (minhash_t) vector_prod >> (w - M);
	}
};
typedef std::vector<rand_hash_function_t> VectorHashFunctions;

struct rand_range_generator_t {
	int rand_in_range(int n) {
		int r, rand_max = RAND_MAX - (RAND_MAX % n);
		while ((r = rand()) >= rand_max);
		return r / (rand_max / n);
	}
};

struct kmer_hasher_t {

	minhash_t encrypt_crypto_32bit(uint32 packed_kmer) const {
		minhash_t cipher = 0;

		return cipher;
	}

	minhash_t encrypt_base_seq(const char* seq, const seq_t seq_len) const {
		return CityHash32(seq, seq_len);
	}
};

inline int irand(int n) {
	int r, rand_max = RAND_MAX - (RAND_MAX % n);
	// reroll until r falls in a range that can be evenly
	// distributed in n bins.  Unless n is comparable to
	// to RAND_MAX, it's not *that* important really. /
	while ((r = rand()) >= rand_max);
	return r / (rand_max / n);
}

// permutes a 64-bit integer (2bit shuffle)
// perm is a random permutation of integers 0..32
inline void perm64(uint64* m, const int* perm) {
	uint64 n = *m;
	uint64_t p = 0;
	//for(int i = 0; i < 64; i++) {
	//	int idx = perm[i];
	//	p |= (((*n >> (idx)) & 1) << (i));
		//p |= (((*n >> (2*idx+1)) & 1) << (2*i+1));
	//}
	p |= (((n >> (perm[0])) & 1) << (0));
	p |= (((n >> (perm[1])) & 1) << (1));
	p |= (((n >> (perm[2])) & 1) << (2));
	p |= (((n >> (perm[3])) & 1) << (3));
	p |= (((n >> (perm[4])) & 1) << (4));
	p |= (((n >> (perm[5])) & 1) << (5));
	p |= (((n >> (perm[6])) & 1) << (6));
	p |= (((n >> (perm[7])) & 1) << (7));
	p |= (((n >> (perm[8])) & 1) << (8));
	p |= (((n >> (perm[9])) & 1) << (9));
	p |= (((n >> (perm[10])) & 1) << (10));
	p |= (((n >> (perm[11])) & 1) << (11));
	p |= (((n >> (perm[12])) & 1) << (12));
	p |= (((n >> (perm[13])) & 1) << (13));
	p |= (((n >> (perm[14])) & 1) << (14));
	p |= (((n >> (perm[15])) & 1) << (15));
	p |= (((n >> (perm[16])) & 1) << (16));
	p |= (((n >> (perm[17])) & 1) << (17));
	p |= (((n >> (perm[18])) & 1) << (18));
	p |= (((n >> (perm[19])) & 1) << (19));
	p |= (((n >> (perm[20])) & 1) << (20));
	p |= (((n >> (perm[21])) & 1) << (21));
	p |= (((n >> (perm[22])) & 1) << (22));
	p |= (((n >> (perm[23])) & 1) << (23));
	p |= (((n >> (perm[24])) & 1) << (24));
	p |= (((n >> (perm[25])) & 1) << (25));
	p |= (((n >> (perm[26])) & 1) << (26));
	p |= (((n >> (perm[27])) & 1) << (27));
	p |= (((n >> (perm[28])) & 1) << (28));
	p |= (((n >> (perm[29])) & 1) << (29));
	p |= (((n >> (perm[30])) & 1) << (30));
	p |= (((n >> (perm[31])) & 1) << (31));
	p |= (((n >> (perm[32])) & 1) << (32));
	p |= (((n >> (perm[33])) & 1) << (33));
	p |= (((n >> (perm[34])) & 1) << (34));
	p |= (((n >> (perm[35])) & 1) << (35));
	p |= (((n >> (perm[36])) & 1) << (36));
	p |= (((n >> (perm[37])) & 1) << (37));
	p |= (((n >> (perm[38])) & 1) << (38));
	p |= (((n >> (perm[39])) & 1) << (39));
	p |= (((n >> (perm[40])) & 1) << (40));
	*m = p;
}

inline void shuffle(int *perm) {
	int tmp;
	int len = 32;
	while(len) {
		int j = irand(len);
		if (j != len - 1) {
			tmp = perm[j];
			perm[j] = perm[len-1];
			perm[len-1] = tmp;
		}
		len--;
	}
}

#endif /*HASH_H_*/
