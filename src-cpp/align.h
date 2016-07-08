#ifndef ALIGN_H_
#define ALIGN_H_

#include "io.h"
#include "index.h"
#include "sam.h"

void balaur_main(const char* fastaName,ref_t& ref, reads_t& reads, precomp_contig_io_t& contig_io, sam_writer_t& sam_io);
void eval(reads_t& reads, const ref_t& ref) ;

#endif /*ALIGN_H_*/
