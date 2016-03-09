#ifndef CONTIGS_H
#define CONTIGS_H

#pragma once
#include "types.h"
#include "index.h"

#define N_TABLES_MAX 1024
#define CONTIG_PADDING 50
#define MAX_BUCKET_SIZE 1000

void assemble_candidate_contigs(const ref_t& ref, reads_t& reads);

#endif