#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "io.h"
#include "cluster.h"

void cluster_reads(reads_t* reads, clusters_t** out) {
	
	clusters_t* clusters = (clusters_t*) malloc(sizeof(clusters_t));
	clusters->num_clusters = 0;
	clusters->alloc_num_clusters = INIT_NUM_CLUSTERS;
	clusters->clusters = (cluster_t*) malloc(clusters->alloc_num_clusters * sizeof(cluster_t));
	
	for(int i = 0; i < reads->count; i++) {
		read_t r = reads->reads[i];	
		int mapped = 0;
		for(int j = 0; j < clusters->num_clusters; j++) {
			cluster_t* cluster = &clusters->clusters[j];			
			if(r.simhash == cluster->simhash) {
				if(cluster->alloc_size == cluster->size) {
					cluster->alloc_size <<= 1;
					cluster->reads = (read_t**) realloc(cluster->reads, cluster->alloc_size * sizeof(read_t*));
				}
				cluster->reads[cluster->size] = &r;
				cluster->size++;
				mapped = 1;
				break;
			} 
		}
		if(!mapped) {
			// assign it to a new cluster
			if(clusters->alloc_num_clusters == clusters->num_clusters) {
				clusters->alloc_num_clusters <<= 1;
				clusters->clusters = (cluster_t*) realloc(clusters->clusters, clusters->alloc_num_clusters * sizeof(cluster_t));
			}
			cluster_t* new_cluster = &clusters->clusters[clusters->num_clusters];
			new_cluster->reads =  (read_t**) malloc(INIT_CLUSTER_SIZE * sizeof(read_t*));
			new_cluster->alloc_size = INIT_CLUSTER_SIZE;
			new_cluster->simhash = r.simhash;
			new_cluster->reads[0] = &r;
			new_cluster->size = 1;
			clusters->num_clusters++;
		}
	}
	*out = clusters;
}

