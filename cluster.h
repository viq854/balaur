#ifndef CLUSTER_H_
#define CLUSTER_H_
#include "io.h"
#include "index.h"

#define REF_MATCHES	1
#define INIT_CLUSTER_SIZE 100
#define INIT_NUM_CLUSTERS 10000

// cluster of similar reads
typedef struct {
	simhash_t simhash; // simhash fp of the representative read 
	int size; // number of reads in the cluster
	int alloc_size;
	read_t** reads; // members of the cluster
	
	seq_t* ref_matches;
	int num_matches;
	int alloc_matches;
	char acc; // DEBUG: whether read matched accuretely
	int best_hamd;
	int best_pos;
} cluster_t;

// clusters of similar reads
typedef struct {
	int num_clusters; // number of clusters
	int alloc_num_clusters;
	cluster_t* clusters;
} clusters_t;

void sort_windows_simhash(ref_t* ref);
void sort_reads_simhash(reads_t* reads);
seq_t get_num_distinct(ref_t* ref);
void cluster_sorted_reads(reads_t* reads, clusters_t** out);
void cluster_reads(reads_t* reads, clusters_t** out);

int collapse_clusters(clusters_t* clusters, index_params_t* params);

#endif /*CLUSTER_H_*/
