#include "contigs.h"
#include <bitset>
#include <string.h>

// Priority heap support (used to process reference index buckets)
struct heap_entry_t {
	seq_t pos;
	uint32_t len;
	uint16_t tid;
	uint16_t next_idx;
};

struct heap_ops {	
	static void heap_sort(heap_entry_t* heap, int n) {
		heap_entry_t tmp;
		int i, j;
		for(j = 1; j < n; j++) {
			tmp = heap[j];
			for(i = j - 1; (i >= 0) && (heap[i].pos > tmp.pos); i--) {
				heap[i+1] = heap[i];
			}
			heap[i+1] = tmp;
		}
	}
	
	static inline void heap_set_min(heap_entry_t* heap, int n) {
		int min_idx = 0;
		for(int i = 1; i < n; i++) {
			if(heap[i].pos < heap[min_idx].pos) {
				min_idx = i;
			}
		}
		heap_entry_t tmp = heap[0];
		heap[0] = heap[min_idx];
		heap[min_idx] = tmp;
	}
	
	static inline void heap_update(heap_entry_t* heap, uint32 n) {
		uint32 i = 0;
		uint32 k = i;
		heap_entry_t tmp = heap[i];
		while((k = (k << 1) + 1) < n) {
			if(k != (n - 1) && (heap[k].pos >= heap[k+1].pos)) ++k;
			if(heap[k].pos >= tmp.pos) break;
			heap[i] = heap[k];
			i = k;
		}
		heap[i] = tmp;
	}
	
	static void swap(heap_entry_t*x, heap_entry_t*y) {
		heap_entry_t temp = *x;  *x = *y;  *y = temp;
	}
	
	static void sift_down(heap_entry_t* heap, int n, int i) {
		int l = 2*i + 1;
		int r = 2*i + 2;
		int min = i;
		if (l < n && heap[l].pos < heap[i].pos) {
			min = l;
		}
		if (r < n && heap[r].pos < heap[min].pos) {
			min = r;
		}
		if (min != i) {
			swap(&heap[i], &heap[min]);
			sift_down(heap, n, min);
		}
	}
	
	static void heap_create(heap_entry_t* heap, int n) {
		int i = (n-1)/2;
		while(i >= 0) {
			sift_down(heap, n, i);
			i--;
		}
	}

	static inline void heap_update_memmove(heap_entry_t* heap, uint32 n) {
		if(n <= 1) return;
		if(heap[1].pos >= heap[0].pos) return; // nothing to shift

		heap_entry_t tmp = heap[0];
		if(heap[n-1].pos <= heap[0].pos) { // shift all
			memmove(heap, heap+1, (n-1)*sizeof(heap_entry_t));
			heap[n-1] = tmp;
		} else { // binary search: find the first element > tmp
			int i, j, k;
			i = 1;
			j = n;
			while(i != j) {
				k = (i + j)/2;
				if (heap[k].pos < tmp.pos) {
					i = k + 1;
				} else {
					j = k;
				}
			}
			memmove(heap, heap+1, (i)*sizeof(heap_entry_t));
			heap[i-1] = tmp;
		}
	}
};

// find the next smallest hit position
int get_next_contig(const ref_t& ref, const std::vector<std::pair<uint64, minhash_t> >& ref_bucket_matches_by_table, uint32 t, heap_entry_t* entry) {
	const minhash_t read_proj_hash = ref_bucket_matches_by_table[t].second;
	const uint64 bid = ref_bucket_matches_by_table[t].first;
	if(bid == ref.index.bucket_offsets.size()) { // table ignored
		return 0;
	}
	const uint64 bucket_data_offset = ref.index.bucket_offsets[bid];
	const uint64 bucket_data_size = ref.index.bucket_offsets[bid+1] - bucket_data_offset;
	// get the next entry in the bucket that matches the read projection hash value
	bool first = entry->next_idx == 0;
	if(first) {
		loc_t l;
		l.hash = read_proj_hash;
		l.pos = 0;
		l.len = 0;
		std::vector<loc_t>::const_iterator range_start = std::lower_bound(ref.index.buckets_data.begin() + bucket_data_offset, ref.index.buckets_data.begin() + bucket_data_offset + bucket_data_size, l, comp_loc()); 	
		entry->next_idx = std::distance(ref.index.buckets_data.begin() + bucket_data_offset, range_start);
	}
	if(ref.index.buckets_data[bucket_data_offset + entry->next_idx].hash == read_proj_hash) {
		entry->pos = ref.index.buckets_data[bucket_data_offset + entry->next_idx].pos;
		entry->len = ref.index.buckets_data[bucket_data_offset + entry->next_idx].len;
		entry->tid = t;
		entry->next_idx++;
		return 1;
	}
	return 0;
}

void process_contig(ref_match_t contig, read_t* r) {
	// filters
	if(contig.len > params->max_matched_contig_len) return;
	if(contig.n_diff_bucket_hits < (int) params->min_n_hits) return;
	if(contig.n_diff_bucket_hits < (int) (r->best_n_bucket_hits - params->dist_best_hit)) return;

	// passed filters
	if(contig.n_diff_bucket_hits > r->best_n_bucket_hits) { // if more hits than best so far
		r->best_n_bucket_hits = contig.n_diff_bucket_hits;
	}
	contig.pos = (contig.pos >= CONTIG_PADDING) ? contig.pos - CONTIG_PADDING : 0;
	contig.len += 2*CONTIG_PADDING + r->len;
	r->ref_matches.push_back(contig);
	r->n_proc_contigs++;
}

inline void test_and_set_valid_contig(ref_match_t& ref_contig, const int n_proc_contigs, const int max_bucket_hits) {
	ref_contig.valid = true;
	if ((ref_contig.n_diff_bucket_hits < 2 && n_proc_contigs > params->proc_contigs_thr) 
		|| (ref_contig.n_diff_bucket_hits < (int) (max_bucket_hits - params->dist_best_hit)))  {
		ref_contig.valid = false;
	}
}

// output matches (ordered by the number of projections matched)
void find_candidate_contigs(const ref_t& ref, read_t* r, const bool rc) {
	r->ref_matches.reserve(10);

	// priority heap of matched positions
	heap_entry_t heap[params->n_tables];
	int heap_size = 0;
	// push the first entries in each sorted bucket onto the heap
	for(uint32 t = 0; t < params->n_tables; t++) { // for each table
		heap[heap_size].next_idx = 0;
		if(get_next_contig(ref, (rc ? r->ref_bucket_matches_by_table_rc : r->ref_bucket_matches_by_table_f), t, &heap[heap_size]) > 0) {
			heap_size++;
		}
	}
	if(heap_size == 0) return; // all the matched buckets are empty
	heap_ops::heap_sort(heap, heap_size); // build heap

	int n_diff_table_hits = 0;
	int len = 0;
	seq_t last_pos = -1;
	std::bitset<N_TABLES_MAX> occ;
	while(heap_size > 0) {
		heap_entry_t e = heap[0]; // get min
		seq_t e_last_pos = e.pos + e.len - 1;
		if(last_pos == (seq_t) -1 || (e.pos <= last_pos)) {
			if(last_pos == (seq_t) -1) { // first contig
				len = e.len - 1;
				last_pos = e_last_pos;
			} else if(last_pos < e_last_pos) { // extending contig
				len += e_last_pos - last_pos;
				last_pos = e_last_pos;
			}
			if(!occ.test(e.tid)) {
				n_diff_table_hits++;
			}
			occ.set(e.tid);
		} else {
			// found a boundary, store/handle last contig
			ref_match_t contig(last_pos - len + 1, len, rc, n_diff_table_hits);
			process_contig(contig, r);

			// start a new contig
			n_diff_table_hits = 1;
			len = e.len - 1;
			last_pos = e_last_pos;
			occ.reset();
			occ.set(e.tid);
		}
		// push the next match from this bucket
		if(get_next_contig(ref, (rc ? r->ref_bucket_matches_by_table_rc : r->ref_bucket_matches_by_table_f), e.tid, &heap[0]) > 0) {
			heap_ops::heap_update(heap, heap_size);
		} else { // no more entries in this bucket
			heap[0] = heap[heap_size - 1];
			heap_size--;
			heap_ops::heap_update(heap, heap_size);
		}
	}
	// add the last position
	if(last_pos != (seq_t) -1) {
		ref_match_t contig(last_pos - len + 1, len, rc, n_diff_table_hits);
		process_contig(contig, r);
	}
}

///// project and merge the resulting buckets
void assemble_candidate_contigs(const ref_t& ref, reads_t& reads) {
	//#pragma omp parallel for	
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(r->valid_minhash_f) {
			r->ref_bucket_matches_by_table_f.resize(params->n_tables);
			for(uint32 t = 0; t < params->n_tables; t++) {
				const minhash_t proj_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes_f, params->sketch_proj_indices, t*params->sketch_proj_len);
				const uint64_t bid = t*params->n_buckets + params->sketch_proj_hash_func.bucket_hash(proj_hash);
				const uint32 bucket_size = ref.index.bucket_offsets[bid + 1] - ref.index.bucket_offsets[bid];
				if(bucket_size > MAX_BUCKET_SIZE) {
					r->ref_bucket_matches_by_table_f[t] = std::pair<uint64, minhash_t>(ref.index.bucket_offsets.size(), 0);
					continue;
				}
				r->ref_bucket_matches_by_table_f[t] = std::pair<uint64, minhash_t>(bid, proj_hash);
				//_mm_prefetch((const void *)&ref.index.buckets_data[ref.index.bucket_offsets[bid]],_MM_HINT_T0);
			}
			find_candidate_contigs(ref, r, false);
			r->n_match_f = r->ref_matches.size();
		}
		if(r->valid_minhash_rc) {
			r->ref_bucket_matches_by_table_rc.resize(params->n_tables);
			for(uint32 t = 0; t < params->n_tables; t++) {
				const minhash_t proj_hash = params->sketch_proj_hash_func.apply_vector(r->minhashes_rc, params->sketch_proj_indices, t*params->sketch_proj_len);
				const uint64_t bid = t*params->n_buckets + params->sketch_proj_hash_func.bucket_hash(proj_hash);
				uint32 bucket_size = ref.index.bucket_offsets[bid + 1] - ref.index.bucket_offsets[bid];
				if(bucket_size > MAX_BUCKET_SIZE) {
					r->ref_bucket_matches_by_table_rc[t] = std::pair<uint64, minhash_t>(ref.index.bucket_offsets.size(), 0);
					continue;
				}
				r->any_bucket_hits = true;
				r->ref_bucket_matches_by_table_rc[t] = std::pair<uint64, minhash_t>(bid, proj_hash);
				//_mm_prefetch((const char *)&ref.index.buckets_data[ref.index.bucket_offsets[bid]],_MM_HINT_T0);
			}
			find_candidate_contigs(ref, r, true);
		}
	}
}

void filter_candidate_contigs(reads_t& reads) {
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->is_valid()) continue;
		bool first_rc = true;
		for(size_t c = 0; c < r->ref_matches.size(); c++) {
			if(r->ref_matches[c].rc && first_rc) {
				r->n_match_f = c;
				first_rc = false;
			} 
			test_and_set_valid_contig(r->ref_matches[c], r->n_proc_contigs, r->best_n_bucket_hits);
		}
		if(first_rc) {
			r->n_match_f = r->ref_matches.size();
		}
	}
}
