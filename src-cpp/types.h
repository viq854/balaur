#ifndef TYPES_H_
#define TYPES_H_

#pragma once

#include <vector>
#include <map>
#include <string>

typedef unsigned int uint32;
typedef unsigned long long int uint64;
typedef unsigned char uint8;

typedef uint64 hash_t;
typedef uint32 minhash_t;
typedef uint32 seq_t;
typedef uint16_t len_t;


typedef uint64 kmer_cipher_t;
typedef seq_t pos_cipher_t;
typedef seq_t pos_offset_cipher_t;

#define MAX_LOC_LEN (1<<16)

struct loc_t {
	seq_t pos;
	len_t len;
	uint32_t hash;
};

struct comp_loc
{
    bool operator()(const loc_t& a, const loc_t& b) const {
    	return a.pos < b.pos;
    }
};

typedef std::vector<uint32> VectorU32;
typedef std::vector<uint8> VectorU8;
typedef std::vector<bool> VectorBool;
typedef std::vector<hash_t> VectorHash;
typedef std::vector<minhash_t> VectorMinHash;

#if(USE_TBB)
#include <tbb/tbb.h>
#include "tbb/scalable_allocator.h"
typedef std::vector<loc_t, tbb::scalable_allocator<loc_t>> VectorSeqPos;
#else
typedef std::vector<loc_t> VectorSeqPos;
#endif

#if(USE_MARISA)
#include <marisa.h>
typedef marisa::Trie MarisaTrie;
#else
struct none_t {};
typedef none_t MarisaTrie;
#endif

#endif /*TYPES_H_*/
