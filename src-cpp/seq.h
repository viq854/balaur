#ifndef SEQ_H
#define SEQ_H

#pragma once
#include "types.h"

// compression
#define CHARS_PER_SHORT 8   // number of chars in 16 bits
#define CHARS_PER_WORD 	16	// number of chars in 32 bits
#define BITS_PER_CHAR 	2
#define BITS_IN_SHORT 	16
#define BITS_IN_WORD 	32
#define BITS_IN_LWORD 	64
#define BASE_IGNORE		4
#define KMER_HIST_SIZE16 (1ULL << 16) //65536
#define KMER_HIST_SIZE32 (1ULL << 32)

int pack_16(const char *seq, const int length, uint16_t *ret);
int pack_32(const char *seq, const int length, uint32_t *ret); 
int pack_64(const char *seq, const int length, uint64 *ret);
void unpack_32(uint32 w, unsigned char *seq, const uint32 length);

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
	bool allowN;

	void init(const std::string& seq, int k) {
		s = seq;
		kmer_len = k;
		pos = 0;
		first = true;
		allowN = false;
	}

	// compression
	inline void pack_init() {
 		kmer.valid = true;
		kmer.packed = 0;
		for (int i = 0; i < kmer_len; i++) {
			uint8 c = s[pos+i];
			if(c == BASE_IGNORE && !allowN) {
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
		if(c == BASE_IGNORE && !allowN) {
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

static inline int get_n_sampled_kmers(const int len, const int k, const int sampling_ratio) {
	return get_n_kmers(len, k)/sampling_ratio;
}

static void find_repeats(const std::string& seq, const int k, std::vector<bool> repeat_mask) {
	const int n_kmers = get_n_kmers(seq.size(), k);
	std::vector<std::pair<packed_kmer_t, int>> kmers(n_kmers);
	repeat_mask.resize(n_kmers);
	kmer_parser_t seq_parser;
	seq_parser.init(seq, k);
	seq_parser.allowN = true;
	kmer_t kmer;
	for(int i = 0; i < n_kmers; i++) {
		seq_parser.get_next_kmer(kmer);
		kmers[i] = std::make_pair(kmer.packed, i);
	}
	std::sort(kmers.begin(), kmers.end());
	
	for(int i = 1; i < n_kmers; i++) {
		if(kmers[i].first == kmers[i-1].first) {
			repeat_mask[kmers[i].second] = true;
			repeat_mask[kmers[i-1].second] = true;
		}
	}
}

#endif