#include "index.h"

void generate_sha1_ciphers(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_len, const std::vector<bool>& repeat_mask, bool rev_mask);
void generate_vanilla_ciphers(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_len);
void apply_keys(kmer_cipher_t* ciphers, const int n_ciphers, const uint64 key1, const uint64 key2);
void mask_repeats(kmer_cipher_t* ciphers, const int n_ciphers);
void  lookup_sha1_ciphers(kmer_cipher_t* ciphers, const seq_t offset, const seq_t len, const std::vector<kmer_cipher_t>& precomp_ref_hashes, const std::vector<uint16_t>& repeat_info);
void  lookup_vanilla_ciphers(kmer_cipher_t* ciphers, const seq_t offset, const seq_t len, const std::vector<kmer_cipher_t>& precomp_ref_hashes);