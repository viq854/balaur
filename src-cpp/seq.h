#ifndef SEQ_H
#define SEQ_H

#pragma once
#include "types.h"
#include "io.h"

#define CHAR_MASK 3ULL
typedef uint32 packed_kmer_t;
struct kmer_t {
	packed_kmer_t packed;
	bool valid; // does not contain any ambiguous bases
};

struct kmer_parser_t {
	std::string s; // sequence to parse
	int kmer_len;
	seq_t pos; // position in the sequence
	kmer_t kmer;
	bool first;

	void init(const std::string& seq, int k) {
		s = seq;
		kmer_len = k;
		pos = 0;
		first = true;
	}

	// compression
	inline void pack_init() {
 		kmer.valid = true;
		kmer.packed = 0;
		for (int i = 0; i < kmer_len; i++) {
			uint8 c = s[pos+i];
			if(c == BASE_IGNORE) {
				kmer.valid = false;
				break;
			}
			kmer.packed |= ((c & CHAR_MASK) << (BITS_IN_WORD - (i+1) * BITS_PER_CHAR));
		}
		if(kmer.valid) {
			pos += kmer_len;
		} else {
			pos++;
		}
	}

	void pack_roll() {
		uint8 c = s[pos];
		if(c == BASE_IGNORE) {
			kmer.valid = false;
			pos++; // skip all the kmers containing this char
			return;
		}
		kmer.packed <<= BITS_PER_CHAR;
		kmer.packed |= ((c & CHAR_MASK) << (BITS_IN_WORD - kmer_len * BITS_PER_CHAR));
		pos++;
	}

	bool get_next_kmer(kmer_t& new_kmer) {
		if(pos >= s.size()) return false;
		// process the next char in the sequence
		if (first || !kmer.valid) {
			pack_init();
			if(first) {
				first = false;
			}
		} else {
			pack_roll();
		}
		new_kmer = kmer;
		return true;
	}
};

static inline int get_n_kmers(const int len, const int k) {
	return len - k + 1;
}

#endif