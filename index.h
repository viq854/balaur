#ifndef INDEX_H_
#define INDEX_H_
#include "io.h"

typedef struct {
	int k; // length of the k-mers
	
	float min_freq;
	float max_freq;
	uint64_t min_count;
	uint64_t max_count;
	
	int max_hammd; // maximum hamming distance to
	
	int ref_window_size;
	
	int s; // length of the hash vector
} index_params_t;

void index_ref(char* fastaFname, index_params_t* params, ref_t** refidx);
void index_reads(char* readsFname, ref_t* ref, index_params_t* params, reads_t** ridx);

#endif /*INDEX_H_*/
