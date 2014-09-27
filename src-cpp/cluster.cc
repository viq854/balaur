//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
//#include <assert.h>
//#include <stdint.h>
//#include <limits.h>
//#include <math.h>
//#include <time.h>
//#include <algorithm>
//
//#include "io.h"
//#include "cluster.h"
//#include "hash.h"
//
//// ref window comparator
//int comp_windows(const ref_win_t r1, const ref_win_t r2) {
//	hash_t h1 = r1.simhash;
//	hash_t h2 = r2.simhash;
//	return (h1 > h2) - (h1 < h2);
//}
//
//// read comparator
//int comp_reads(const read_t r1, const read_t r2) {
//	hash_t h1 = r1.simhash;
//	hash_t h2 = r2.simhash;
//	return (h1 > h2) - (h1 < h2);
//}
//
//// cluster comparator
//int comp_clusters(const cluster_t r1, const cluster_t r2) {
//	hash_t h1 = r1.simhash;
//	hash_t h2 = r2.simhash;
//	return (h1 > h2) - (h1 < h2);
//}
//
//// sorts reference windows by their simhash value
//void sort_windows_hash(ref_t& ref) {
//	//qsort(ref->windows, ref->num_windows, sizeof(ref_win_t), comp_windows);
//}
//
//// sorts clusters by their simhash value
//void sort_clusters_hash(VectorClusters& clusters) {
//	std::sort(clusters.begin(), clusters.end(), comp_clusters);
//}
//
//// sorts reads by their simhash value
//void sort_reads_hash(reads_t& reads) {
//	std::sort(reads.reads.begin(), reads.reads.end(), comp_reads);
//}
//
//// finds the number of windows with a different simhash
//seq_t get_num_distinct(ref_t& ref) {
//	//hash_t prev;
//	seq_t num_diff = 0;
////	for(MapPos2Window::iterator it = ref.windows_by_pos.begin(); it != ref.windows_by_pos.end(); ++it) {
////		ref_win_t w = it->second;
////		if((it == ref.windows_by_pos.begin()) || (w.simhash != prev)) {
////			num_diff++;
////			prev = w.simhash;
////		}
////	}
//	return num_diff;
//}
//
//// collapse clusters with similar simhash value
//#define MAX_DIST 65
//int collapse_clusters(VectorClusters& clusters, index_params_t* params) {
//	sort_clusters_hash(clusters);
//
//	// find the Hamming distance between adjacent simhashes (TODO: consider a window)
//	VectorU32 hammd_pairs(clusters.size() - 1);
//	uint32 min_dist = MAX_DIST;
//	uint32 min_idx = 0;
//	for(uint32 i = 0; i < clusters.size() - 1; i++) {
//		hash_t h1 = clusters[i].simhash;
//		hash_t h2 = clusters[i+1].simhash;
//		uint32 dist = hamming_dist(h1, h2);
//		hammd_pairs[i] = dist;
//		if(dist < min_dist) {
//			min_dist = dist;
//			min_idx = i;
//		}
//	}
//
//	int num_collapsed = 0;
//	while(min_dist < params->max_hammd) {
//		// collapse cluster pair and min_idx
//		cluster_t* c1 = &clusters[min_idx];
//		cluster_t* c2 = &clusters[min_idx+1];
//
//		for(uint32 i = 0; i < c2->reads.size(); i++) {
//			c1->reads.push_back(c2->reads[i]);
//		}
//
//		// remove the pair and the neighbors from min search
//		hammd_pairs[min_idx] = MAX_DIST;
//		if(min_idx + 1 < clusters.size()) hammd_pairs[min_idx + 1] = MAX_DIST;
//		if(min_idx - 1 >= 0) hammd_pairs[min_idx - 1] = MAX_DIST;
//
//		// find the next min
//		min_dist = MAX_DIST;
//		for(uint32 i = 0; i < clusters.size() - 1; i++) {
//			if(hammd_pairs[i] < min_dist) {
//				min_dist = hammd_pairs[i];
//				min_idx = i;
//			}
//		}
//		num_collapsed++;
//	}
//	return num_collapsed;
//}
//
//// finds the reads with the same simhash value and assigns them into the same cluster
//void cluster_sorted_reads(reads_t& reads, VectorClusters& clusters) {
//	uint32 prev_cluster_id;
//	for(uint32 i = 0; i < reads.reads.size(); i++) {
//		if((i > 0) && (i % 100000 == 0)) {
//			printf("Clustered %u reads. Total of %zu distinct clusters found.\n", i, clusters.size());
//		}
//		read_t* r = &reads.reads[i];
//		if((i == 0) || (r->simhash != clusters[prev_cluster_id].simhash)) {
//			cluster_t new_cluster;
//			new_cluster.simhash = r->simhash;
//			new_cluster.reads.push_back(r);
//			new_cluster.acc = 0;
//			new_cluster.best_hamd = INT_MAX;
//			clusters.push_back(new_cluster);
//			prev_cluster_id = i;
//		} else {
//			clusters[prev_cluster_id].reads.push_back(r);
//		}
//	}
//}
//
//void cluster_reads(reads_t& reads, VectorClusters& clusters) {
//
//	for(uint32 i = 0; i < reads.reads.size(); i++) {
//		read_t* r = &reads.reads[i];
//		int mapped = 0;
//		for(uint32 j = 0; j < clusters.size(); j++) {
//			cluster_t* cluster = &clusters[j];
//			if(r->simhash == cluster->simhash) {
//				cluster->reads.push_back(r);
//				mapped = 1;
//				break;
//			}
//		}
//		if(!mapped) {
//			// assign it to a new cluster
//
//			cluster_t new_cluster;
//			new_cluster.simhash = r->simhash;
//			new_cluster.reads.push_back(r);
//			clusters.push_back(new_cluster);
//		}
//	}
//}
//
