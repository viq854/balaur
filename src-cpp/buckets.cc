struct heap_entry_t {
	seq_t pos;
	uint32_t len;
	uint16_t tid;
	uint16_t next_idx;
};

void heap_sort(heap_entry_t* heap, int n) {
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

inline void heap_set_min(heap_entry_t* heap, int n) {
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

inline void heap_update(heap_entry_t* heap, uint32 n) {
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

// min-heap
void swap(heap_entry_t*x, heap_entry_t*y) {
	heap_entry_t temp = *x;  *x = *y;  *y = temp;
}

void sift_down(heap_entry_t* heap, int n, int i) {
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

void heap_create(heap_entry_t* heap, int n) {
	int i = (n-1)/2;
	while(i >= 0) {
		sift_down(heap, n, i);
		i--;
	}
}

inline void heap_update_memmove(heap_entry_t* heap, uint32 n) {
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

int get_next_contig(const ref_t& ref, const std::vector<std::pair<uint64, minhash_t> >& ref_bucket_matches_by_table, uint32 t, heap_entry_t* entry) {
	const minhash_t read_proj_hash = ref_bucket_matches_by_table[t].second;
	const uint64 bid = ref_bucket_matches_by_table[t].first;
	if(bid == ref.index.bucket_offsets.size()) { // table ignored
		return 0;
	}
	const uint64 bucket_data_offset = ref.index.bucket_offsets[bid];
	const uint64 bucket_data_size = ref.index.bucket_offsets[bid+1] - bucket_data_offset;

	//const uint64 bucket_data_offset = r->bucket_offsets[t];
	//const uint64 bucket_data_size = r->bucket_offsets[t+1] - bucket_data_offset;

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
	//while(entry->next_idx < bucket_data_size) {
		if(ref.index.buckets_data[bucket_data_offset + entry->next_idx].hash == read_proj_hash) {
			entry->pos = ref.index.buckets_data[bucket_data_offset + entry->next_idx].pos;
			entry->len = ref.index.buckets_data[bucket_data_offset + entry->next_idx].len;
			entry->tid = t;
			entry->next_idx++;
			return 1;
		}
		//if(!first) {
		//	return 0;
		//}
		//entry->next_idx++;
	//}
	return 0;
}