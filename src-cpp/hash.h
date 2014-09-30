#ifndef HASH_H_
#define HASH_H_

#pragma once
#include "types.h"
#include "city.h"
#include "mt64.h"
#include "../rollinghashcpp/cyclichash.h"

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
		return (minhash_t) s >> (w - M);
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

void hashlittle2( 
  const void *key,       /* the key to hash */
  size_t      length,    /* length of the key */
  uint32_t   *pc,        /* IN: primary initval, OUT: primary hash */
  uint32_t   *pb); 

#endif /*HASH_H_*/
