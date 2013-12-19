#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "io.h"
#include "cluster.h"

// read comparator
int comp_reads(const void * r1, const void * r2) {
	simhash_t h1 = ((read_t*) r1)->simhash;
	simhash_t h2 = ((read_t*) r2)->simhash;
	return (h1 > h2) - (h1 < h2);
}

// cluster comparator
int comp_clusters(const void * r1, const void * r2) {
	simhash_t h1 = ((read_t*) r1)->simhash;
	simhash_t h2 = ((read_t*) r2)->simhash;
	return (h1 > h2) - (h1 < h2);
}

// sorts reads by their simhash value
void sort_reads_simhash(reads_t* reads) {
	qsort(reads->reads, reads->count, sizeof(read_t), comp_reads);
}

// sorts clusters by their simhash value
void sort_clusters_simhash(clusters_t* clusters) {
	qsort(clusters->clusters, clusters->num_clusters, sizeof(cluster_t), comp_clusters);
}

// collapse clusters with similar simhash value
void collapse_clusters(clusters_t* clusters) {
	
}

// finds the reads with the same simhash value and assigns them into the same cluster
void cluster_sorted_reads(reads_t* reads, clusters_t** out) {
	clusters_t* clusters = (clusters_t*) malloc(sizeof(clusters_t));
	clusters->num_clusters = 0;
	clusters->alloc_num_clusters = INIT_NUM_CLUSTERS;
	clusters->clusters = (cluster_t*) malloc(clusters->alloc_num_clusters * sizeof(cluster_t));
	
	cluster_t* prev_cluster;
	for(int i = 0; i < reads->count; i++) {
		if((i > 0) && (i % 50000 == 0)) {
			printf("Clustered %d reads. Total of %d distinct clusters found.\n", i, clusters->num_clusters);
		}
		read_t* r = &reads->reads[i];	
		if((i == 0) || (r->simhash != prev_cluster->simhash)){
			// assign it to a new cluster
			if(clusters->alloc_num_clusters == clusters->num_clusters) {
				clusters->alloc_num_clusters <<= 1;
				clusters->clusters = (cluster_t*) realloc(clusters->clusters, clusters->alloc_num_clusters * sizeof(cluster_t));
			}
			cluster_t* new_cluster = &clusters->clusters[clusters->num_clusters];
			new_cluster->reads =  (read_t**) malloc(INIT_CLUSTER_SIZE * sizeof(read_t*));
			new_cluster->alloc_size = INIT_CLUSTER_SIZE;
			new_cluster->simhash = r->simhash;
			new_cluster->reads[0] = r;
			new_cluster->size = 1;
			clusters->num_clusters++;
			prev_cluster = new_cluster;
		} else {
			if(prev_cluster->alloc_size == prev_cluster->size) {
				prev_cluster->alloc_size <<= 1;
				prev_cluster->reads = (read_t**) realloc(prev_cluster->reads, prev_cluster->alloc_size * sizeof(read_t*));
			}
			prev_cluster->reads[prev_cluster->size] = r;
			prev_cluster->size++;
		}
	}
	*out = clusters;
}

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

