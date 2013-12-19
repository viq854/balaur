#ifndef INDEX_H_
#define INDEX_H_

typedef struct {
	int k; // length of the k-mers
	float min_freq;
	float max_freq;
	uint64_t min_count;
	uint64_t max_count;
	
	int max_hammd; // maximum hamming distance to
	
	int s; // length of the hash vector
} index_params_t;

void index_reads(char* readsFname, index_params_t* params);

#endif /*INDEX_H_*/
