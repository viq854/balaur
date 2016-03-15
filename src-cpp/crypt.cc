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
				ciphers[i]  = genrand64_int64();
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

// contig hashing
static int bin_shuffle[20] = {0,2,4,6,8,10,12,14,16,18,1,3,5,7,9,11,13,15,17,19};
void lookup_sha1_ciphers(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_offset, const seq_t seq_len, const ref_t& ref) {
	const int n_kmers = get_n_kmers(seq_len, params->k2);
	const int n_sampled_kmers = get_n_sampled_kmers(seq_len, params->k2, params->sampling_intv);
	for(int i = 0; i < n_sampled_kmers; i++) {
		const seq_t p = i*params->sampling_intv;
		ciphers[i] = ref.precomputed_kmer2_hashes[seq_offset + p];
	}
	
	for(int i = 0; i < n_sampled_kmers; i++) {
		const seq_t p = i*params->sampling_intv;
		const uint16_t r = ref.precomputed_neighbor_repeats[seq_offset + p];
		if(r == 0) continue; // unique kmer
		const seq_t next_occ =  p + r;
		if(next_occ >= n_kmers) continue; // repeat is outside the contig
		ciphers[i] = genrand64_int64();
		if(next_occ % params->sampling_intv == 0) {
			ciphers[next_occ/params->sampling_intv] = genrand64_int64();
		}
	}
}

//static int n_prob_reads = 0;
//void generate_voting_kmer_ciphers_ref(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_offset, const seq_t seq_len,
//		const uint64 key1, const uint64 key2, const ref_t& ref) {

	
	//memcpy(&ciphers[0], &ref.precomputed_kmer2_hashes[seq_offset], n_kmers*sizeof(kmer_cipher_t));
	//ciphers[i] = CityHash64(&seq[seq_offset + p], params->k2);
//#if(!VANILLA)

//	const int n_sampled_kmers = get_n_sampled_kmers(seq_len, params->k2, params->sampling_intv);
//	for(int i = 0; i < n_sampled_kmers; i++) {
 //               seq_t p = i*params->sampling_intv;
//		ciphers[i] = ref.precomputed_kmer2_hashes[seq_offset + p];
//	}

//	for(int i = 0; i < n_sampled_kmers; i++) {
//		seq_t p = i*params->sampling_intv;
//		uint16_t r = ref.precomputed_neighbor_repeats[seq_offset + p];
//		if(ciphers[i] != 0 && (r == 0 || r >= (n_kmers-p))) {
//			ciphers[i] = (ciphers[i] ^ key1)*key2;
//		} else {
//			ciphers[i] = genrand64_int64();
  //                      if(r > 0 && r < (n_kmers-p) && (p+r) % params->sampling_intv == 0) ciphers[(p+r)/params->sampling_intv] = 0;
//		}
//	}
	/*const int n_bins = ceil(((float)n_kmers)/params->k2);
	int idx = 0;
	for(int i = 0; i < n_bins; i++) {
		seq_t offset = i*params->k2;
		int n_unique = 0;
		int bin_size = params->k2;
		if(i == n_bins - 1) bin_size = n_kmers - params->k2*i; 
		const int n_kept = bin_size/params->sampling_intv;
		for(int j = 0; j < bin_size; j++) {
			if(n_unique == n_kept) break;
			int k = bin_shuffle[j];
			if(k >= bin_size) continue;
			seq_t pos = offset + k;
			uint16_t r = ref.precomputed_neighbor_repeats[pos + seq_offset];
                	if(ciphers[pos] != 0 && (r == 0 || r >= (n_kmers-i))) {
				ciphers[idx] = (ciphers[pos] ^ key1)*key2; 		
				idx++;
				n_unique++;
			} else {
				if(r > 0 && r < (n_kmers-i)) ciphers[pos + r] = 0;
				//ciphers[idx] = genrand64_int64();//(ciphers[pos] ^ key1)*key2;
                                //idx++;
				//n_unique++;
			}
		}
		for(int x = n_unique; x < n_kept; x++) {
			ciphers[idx] = genrand64_int64();//(ciphers[pos] ^ key1)*key2;
                        idx++;
               	}
		///if(n_unique < n_kept) {
		//	n_prob_reads++;
		//	break;
		}
	}*/
//#endif
//}
