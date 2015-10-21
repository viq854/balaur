#ifndef ALIGN_H_
#define ALIGN_H_

#include "io.h"
#include "index.h"

void align_reads_lsh(ref_t& ref, reads_t& reads, const index_params_t* params);
void align_reads_minhash(ref_t& ref, reads_t& reads, const index_params_t* params);
void align_reads_sampling(ref_t& ref, reads_t& reads, const index_params_t* params);
void balaur_main(const char* fastaName,ref_t& ref, reads_t& reads, const index_params_t* params);

#endif /*ALIGN_H_*/
