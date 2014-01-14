#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "io.h"
#include "cluster.h"
#include "hash.h"

// ref window comparator
int comp_windows(const void * r1, const void * r2) {
	simhash_t h1 = ((ref_win_t*) r1)->simhash;
	simhash_t h2 = ((ref_win_t*) r2)->simhash;
	return (h1 > h2) - (h1 < h2);
}

// read comparator
int comp_reads(const void * r1, const void * r2) {
	simhash_t h1 = ((read_t*) r1)->simhash;
	simhash_t h2 = ((read_t*) r2)->simhash;
	return (h1 > h2) - (h1 < h2);
}

// cluster comparator
int comp_clusters(const void * r1, const void * r2) {
	simhash_t h1 = ((cluster_t*) r1)->simhash;
	simhash_t h2 = ((cluster_t*) r2)->simhash;
	return (h1 > h2) - (h1 < h2);
}

// sorts reference windows by their simhash value
void sort_windows_simhash(ref_t* ref) {
	qsort(ref->windows, ref->num_windows, sizeof(ref_win_t), comp_windows);
}

// sorts clusters by their simhash value
void sort_clusters_simhash(clusters_t* clusters) {
	qsort(clusters->clusters, clusters->num_clusters, sizeof(cluster_t), comp_clusters);
}

// sorts reads by their simhash value
void sort_reads_simhash(reads_t* reads) {
	qsort(reads->reads, reads->count, sizeof(read_t), comp_reads);
}

// finds the number of windows with a different simhash
seq_t get_num_distinct(ref_t* ref) {
	simhash_t prev;
	seq_t num_diff = 0;
	for(seq_t i = 0; i < ref->num_windows; i++) {
		if((i == 0) || (ref->windows[i].simhash != prev)){
			num_diff++;
			prev = ref->windows[i].simhash;
		} 
	}
	return num_diff;
}

void add_read_to_cluster(cluster_t* cluster, read_t* r) {
	if(cluster->alloc_size == cluster->size) {
		cluster->alloc_size <<= 1;
		cluster->reads = (read_t**) realloc(cluster->reads, cluster->alloc_size * sizeof(read_t*));
	}
	cluster->reads[cluster->size] = r;
	cluster->size++;
}

// collapse clusters with similar simhash value
#define MAX_DIST 65
int collapse_clusters(clusters_t* clusters, index_params_t* params) {
	sort_clusters_simhash(clusters);
	
	// find the Hamming distance between adjacent simhashes (TODO: consider a window)
	int* hammd_pairs  = (int*) malloc((clusters->num_clusters - 1) * sizeof(int));
	
	int min_dist = MAX_DIST;
	int min_idx = 0;
	for(int i = 0; i < clusters->num_clusters - 1; i++) {
		simhash_t h1 = clusters->clusters[i].simhash;
		simhash_t h2 = clusters->clusters[i+1].simhash;
		int dist = hamming_dist(h1, h2);
		hammd_pairs[i] = dist;
		if(dist < min_dist) {
			min_dist = dist;
			min_idx = i;
		}
	}
	
	int num_collapsed = 0;
	while(min_dist < params->max_hammd) {
		// collapse cluster pair and min_idx
		cluster_t* c1 = &clusters->clusters[min_idx];
		cluster_t* c2 = &clusters->clusters[min_idx+1];
		
		for(int i = 0; i < c2->size; i++) {
			add_read_to_cluster(c1, c2->reads[i]);
		}
		
		// remove the pair and the neighbors from min search
		hammd_pairs[min_idx] = MAX_DIST;
		if(min_idx + 1 < clusters->num_clusters) hammd_pairs[min_idx + 1] = MAX_DIST;
		if(min_idx - 1 >= 0) hammd_pairs[min_idx - 1] = MAX_DIST;
		
		// find the next min
		min_dist = MAX_DIST;
		for(int i = 0; i < clusters->num_clusters - 1; i++) {
			if(hammd_pairs[i] < min_dist) {
				min_dist = hammd_pairs[i];
				min_idx = i;
			}
		}	
		num_collapsed++;
	}
	return num_collapsed;
}

// finds the reads with the same simhash value and assigns them into the same cluster
void cluster_sorted_reads(reads_t* reads, clusters_t** out) {
	clusters_t* clusters = (clusters_t*) malloc(sizeof(clusters_t));
	clusters->num_clusters = 0;
	clusters->alloc_num_clusters = INIT_NUM_CLUSTERS;
	clusters->clusters = (cluster_t*) malloc(clusters->alloc_num_clusters * sizeof(cluster_t));
	
	cluster_t* prev_cluster;
	for(int i = 0; i < reads->count; i++) {
		if((i > 0) && (i % 100000 == 0)) {
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
			new_cluster->ref_matches = NULL;
			new_cluster->num_matches = 0;
			new_cluster->acc = 0;
			new_cluster->best_hamd = INT_MAX;
			clusters->num_clusters++;
			prev_cluster = new_cluster;
		} else {
			add_read_to_cluster(prev_cluster, r);
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

