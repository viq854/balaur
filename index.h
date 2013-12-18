#ifndef INDEX_H_
#define INDEX_H_

typedef struct {
	int k; // length of the k-mers
	int s; // length of the hash vector
} index_params_t;

void index_reads(char* readsFname, index_params_t* params);

#endif /*INDEX_H_*/
