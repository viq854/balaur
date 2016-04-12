#include <emmintrin.h>
#include <smmintrin.h>
#include <openssl/sha.h>
#include <unordered_map>
#include <unordered_set>
#include "hash.h"
#include "index.h"
#include "seq.h"
#include "crypt.h"

// read hashing
// repeat kmers are masked by default according to the repeat_mask
void generate_sha1_ciphers(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_len, const std::vector<bool>& repeat_mask, bool rev_mask) {
		const int n_kmers = get_n_kmers(seq_len, params->k2);
		uint32_t hash[5];
		for(int i = 0; i < n_kmers; i++) {
			int mask_idx = i;
			if(rev_mask) mask_idx = n_kmers-i-1;
			if(repeat_mask[mask_idx]) {
				ciphers[i]  = genrand64_int64();
			} else {
				sha1_hash(reinterpret_cast<const uint8_t*>(&seq[i]), params->k2, hash);
				ciphers[i] = ((uint64) hash[0] << 32 | hash[1]);
			}
		}
}

// read non-cryptographic hashing
// repeats are allowed
void generate_vanilla_ciphers(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_len) {
	const int n_kmers = get_n_kmers(seq_len, params->k2);
	for(int i = 0; i < n_kmers; i++) {
			ciphers[i] = CityHash64(&seq[i], params->k2);
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

void compute_repeat_mask(const seq_t offset, const int len, const std::vector<uint16_t>& repeat_info, std::vector<bool>& repeat_mask, const int bin_size) {
	for(int i = 0; i < len; i++) {
		const uint16_t r = repeat_info[offset + i]; // distance to closest repeat
		const seq_t next_occ =  i + r;
		if(r == 0 || next_occ >= len) continue; // unique kmer or repeat is outside the contig}
		//const int bin_id = i / bin_size;
		if(!repeat_mask[i]) {
			repeat_mask[i] = true;
			//bin_counts[bin_id]++;
		}
		if(!repeat_mask[next_occ]) {
			repeat_mask[next_occ] = true;
			//bin_counts[bin_id]++;
		}
	}
}

inline bool test_and_set_repeat(const seq_t local_pos, const seq_t offset, const int len, const std::vector<uint16_t>& repeat_info, std::vector<bool>& repeat_mask) {
	const uint16_t r = repeat_info[local_pos + offset]; // distance to closest repeat
	const seq_t next_occ = local_pos + r;
        if(r == 0 || next_occ >= len) return repeat_mask[local_pos];
	repeat_mask[next_occ] = true;
	return true;
}


// strided lookup of precomputed ref kmers (access pattern stored in the shuffle array)
void gather_sha1_ciphers(kmer_cipher_t* ciphers, const std::vector<int>& shuffle, const int shuffle_len, const seq_t offset, const std::vector<kmer_cipher_t>& precomp_ref_hashes) {
	for(int i = 0; i < shuffle_len; i++) {
		ciphers[i] = precomp_ref_hashes[offset + shuffle[i]];
	}
}

// contig hashing
// lookup precomputed sha-1 hashes
// mask repeats
void  lookup_sha1_ciphers(kmer_cipher_t* ciphers, const bool any_repeats, const seq_t offset, const seq_t len, const std::vector<kmer_cipher_t>& precomp_ref_hashes, const std::vector<uint16_t>& repeat_info) {
	const int n_kmers = get_n_kmers(len, params->k2);
	const int n_bins = ceil(((float)n_kmers)/params->bin_size);
	int bin_size = params->bin_size;
	std::vector<bool> repeat_mask(n_kmers);
	
	// uniform sampling
	if(!params->bin_sampling()) {
		compute_repeat_mask(offset, n_kmers, repeat_info, repeat_mask, bin_size);
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
			if(idx >= bin_size || test_and_set_repeat(kmer_offset + idx, offset, n_kmers, repeat_info, repeat_mask)) continue; // index out of range in the last bucket or repeat
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

void  lookup_vanilla_ciphers(kmer_cipher_t* ciphers, const seq_t offset, const seq_t len, const std::vector<kmer_cipher_t>& precomp_ref_hashes) {
	const int n_sampled_kmers = get_n_sampled_kmers(len, params->k2, params->sampling_intv);
	for(int i = 0; i < n_sampled_kmers; i++) {
		seq_t pos = i*params->sampling_intv;
		ciphers[i] = precomp_ref_hashes[offset + pos];
	}
}
