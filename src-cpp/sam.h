#ifndef __SAM_H
#define __SAM_H

#include "index.h"

void store_alns_sam(const reads_t& reads, const ref_t& ref, const index_params_t* params);

#endif
