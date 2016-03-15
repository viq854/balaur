#include "index.h"

void generate_sha1_ciphers(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_len, const std::vector<bool> repeat_mask);
void apply_keys(kmer_cipher_t* ciphers, const int n_ciphers, const uint64 key1, const uint64 key2);
void mask_repeats(kmer_cipher_t* ciphers, const int n_ciphers);
void lookup_sha1_ciphers(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_offset, const seq_t seq_len, const ref_t& ref);