#ifndef HASH_H_
#define HASH_H_

#pragma once
#include "types.h"
#include "city.h"
#include "mt64.h"
struct CyclicHash{ CyclicHash(int x, int y){}; };

#if(USE_SHA1_ASM)
extern "C" 
{
 	void sha1_compress(uint32_t state[5], const uint8_t block[64]);
}
#else
void sha1_compress(uint32_t state[5], const uint8_t block[64]);
#endif

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
