#include <emmintrin.h>
#include <smmintrin.h>
#include <openssl/sha.h>
#include <unordered_map>
#include <unordered_set>
#include "hash.h"
#include "index.h"
#include "seq.h"
#include "crypt.h"

static int max_repeats = 0;
static int max_repeats_contig = 0;
static std::vector<int> n_repeats_v;
static std::vector<int> n_repeats_v_contig;

// read hashing
// repeat kmers are masked by default according to the repeat_mask
void generate_sha1_ciphers(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_len, const std::vector<bool>& repeat_mask, bool rev_mask) {
		const int n_kmers = get_n_kmers(seq_len, params->k2);
		uint32_t hash[5];
		for(int i = 0; i < n_kmers; i++) {
			int mask_idx = i;
			if(rev_mask) mask_idx = n_kmers-i-1;
			if(repeat_mask[mask_idx]) {
				ciphers[i]  = 0; //genrand64_int64();
			} else {
				sha1_hash(reinterpret_cast<const uint8_t*>(&seq[i]), params->k2, hash);
				ciphers[i] = ((uint64) hash[0] << 32 | hash[1]);
			}
		}
}

void apply_keys(kmer_cipher_t* ciphers, const int n_ciphers, const uint64 key1, const uint64 key2) {
	__m128i* c = (__m128i*)ciphers;
	__m128i xor_pad = _mm_set1_epi64((__m64)key1);
	for(int i = 0; i < n_ciphers/2; i++) {
			c[i] = _mm_xor_si128(c[i], xor_pad);
			ciphers[2*i] *= key2;
			ciphers[2*i+1] *= key2;
	}
	for(int i = 2*(n_ciphers/2); i < n_ciphers; i++) {
			ciphers[i] ^= key1;
			ciphers[i] *= key2;
	}
}

// simple repeat masking using a map for lookups
void mask_repeats(kmer_cipher_t* ciphers, const int n_ciphers) {
	std::unordered_map<kmer_cipher_t, int> s;
	std::pair<std::unordered_map<kmer_cipher_t, int>::iterator, bool> r;
	for(int i = 0; i < n_ciphers; i++) {
		r = s.insert(std::make_pair(ciphers[i], i));
		if(!r.second) {
			ciphers[(r.first)->second] = genrand64_int64();
			ciphers[i] = genrand64_int64();
		}
	}
}

void compute_repeat_mask(const seq_t offset, const int len, const std::vector<uint16_t>& repeat_info, std::vector<bool>& repeat_mask) {
	repeat_mask.resize(len);
	for(int i = 0; i < len; i++) {
		const uint16_t r = repeat_info[offset + i]; // distance to closest repeat
		if(r == 0) continue; // unique kmer
		const seq_t next_occ =  i + r;
		if(next_occ >= len) continue; // repeat is outside the contig
		repeat_mask[i] = true;
		repeat_mask[next_occ] = true;
	}
}

// strided lookup of precomputed ref kmers (access pattern stored in the shuffle array)
void gather_sha1_ciphers(kmer_cipher_t* ciphers, const std::vector<int>& shuffle, const int shuffle_len, const seq_t offset, const std::vector<kmer_cipher_t>& precomp_ref_hashes) {
	for(int i = 0; i < shuffle_len; i++) {
		ciphers[i] = precomp_ref_hashes[offset + shuffle[i]];
	}
}

void  lookup_sha1_ciphers(kmer_cipher_t* ciphers, const seq_t offset, const seq_t len, const std::vector<kmer_cipher_t>& precomp_ref_hashes, const std::vector<uint16_t>& repeat_info) {
	const int n_kmers = get_n_kmers(len, params->k2);
	// mark repeats
	std::vector<bool> repeat_mask;
	compute_repeat_mask(offset, n_kmers, repeat_info, repeat_mask);
	
	// uniform sampling
	if(!params->bin_sampling()) {
		if(params->mask_repeat_nbrs) {
			mask_repeat_nbrs(repeat_mask, params->k2);
		}
		const int n_sampled_kmers = get_n_sampled_kmers(len, params->k2, params->sampling_intv);
		for(int i = 0; i < n_sampled_kmers; i++) {
			const int pos = i*params->sampling_intv;
			if(repeat_mask[pos]) {
				ciphers[i] = genrand64_int64();
			} else {
				ciphers[i] = precomp_ref_hashes[offset + pos];
			}
		}
		return;
	}
	
	// bin sampling
	const int n_bins = ceil(((float)n_kmers)/params->bin_size);
	int bin_size = params->bin_size;
	int n_sampled = bin_size/params->sampling_intv;
	std::vector<int> shuffle(n_sampled);
	for(int i = 0; i < n_bins; i++) {
		const int kmer_offset = i*params->bin_size;
		const int cipher_offset = i*n_sampled;
		if(i == n_bins -1)  {
			bin_size = n_kmers - kmer_offset;
			n_sampled = bin_size/params->sampling_intv;
		}
		int n_unique = 0;
		for(int j = 0; j < bin_size; j++) {
			if(n_unique == n_sampled) break; // sampled sufficient unique kmers 
			const int idx = params->bin_shuffle[j];
			if(idx >= bin_size || repeat_mask[kmer_offset +idx]) continue; // index out of range in the last bucket or repeat
			shuffle[n_unique] = idx;
			n_unique++;
		}
		if(n_unique != n_sampled && params->mask_repeat_nbrs) {
			n_unique = 0;
		}

		gather_sha1_ciphers(&ciphers[cipher_offset], shuffle, n_unique, offset + kmer_offset, precomp_ref_hashes); // fill in the sampled hashes
		for(int j = n_unique; j < n_sampled; j++) {
			ciphers[cipher_offset + j] = genrand64_int64();
		}
	}
}