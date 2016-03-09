#include "index.h"

void generate_voting_kmer_ciphers_read(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_len, const uint64 key1, const uint64 key2);
void generate_voting_kmer_ciphers_ref(kmer_cipher_t* ciphers, const char* seq, const seq_t seq_offset, const seq_t seq_len,
		const uint64 key1, const uint64 key2, const ref_t& ref);