#ifndef CLUSTER_H_
#define CLUSTER_H_
#include "io.h"

#define INIT_CLUSTER_SIZE 100
#define INIT_NUM_CLUSTERS 10000

// cluster of similar reads
typedef struct {
	simhash_t simhash; // simhash fp of the representative read 
	int size; // number of reads in the cluster
	int alloc_size;
	read_t** reads; // members of the cluster
} cluster_t;

// clusters of similar reads
typedef struct {
	int num_clusters; // number of clusters
	int alloc_num_clusters;
	cluster_t* clusters;
} clusters_t;

void cluster_reads(reads_t* reads, clusters_t** out);

#endif /*CLUSTER_H_*/
