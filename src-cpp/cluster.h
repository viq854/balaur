#ifndef CLUSTER_H_
#define CLUSTER_H_
#include "io.h"
#include "index.h"

#define REF_MATCHES	1
#define INIT_CLUSTER_SIZE 100
#define INIT_NUM_CLUSTERS 10000

// cluster of similar reads
typedef struct {
	hash_t simhash; 			// simhash fp of the representative read
	VectorPReads reads; 		// members of the cluster
	VectorSeqPos ref_matches;	// match positions in the reference

	char acc; // DEBUG: whether read matched accuretely
	int best_hamd;
	int best_pos;
} cluster_t;
typedef std::vector<cluster_t> VectorClusters;

void sort_windows_hash(ref_t& ref);
void sort_reads_hash(reads_t& reads);
seq_t get_num_distinct(ref_t& ref);
void cluster_sorted_reads(reads_t& reads, VectorClusters& out);
void cluster_reads(reads_t& reads, VectorClusters& out);

int collapse_clusters(VectorClusters& clusters, index_params_t& params);

#endif /*CLUSTER_H_*/
