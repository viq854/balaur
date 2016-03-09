#include "index.h"

static const int VERBOSE = (getenv("VERBOSE") ? atoi(getenv("VERBOSE")) : 0);
void eval(reads_t& reads, const ref_t& ref) {
	printf("////////////// Evaluation //////////////\n");
	// ---- debug -----
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if (VERBOSE > 0 && r->top_aln.score >= 10 &&
				!pos_in_range(r->ref_pos_r, r->top_aln.ref_start, 30) &&
				!pos_in_range(r->ref_pos_l, r->top_aln.ref_start, 30)) {
				printf("WRONG: score %u max-votes: %u second-best-votes: %u true-contig-votes: %u true-bucket-hits: %u max-bucket-hits %u true-pos-l  %u true-pos-r: %u found-pos %u\n",
					r->top_aln.score,
					r->top_aln.inlier_votes, r->second_best_aln.inlier_votes, r->comp_votes_hit, r->true_n_bucket_hits, r->best_n_bucket_hits,
					r->ref_pos_l, r->ref_pos_r, r->top_aln.ref_start);
			//print_read(r);
		}
	}

#if(!SIM_EVAL)
	int q10a = 0;
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->valid_minhash_f && !r->valid_minhash_rc) continue;
		if(r->top_aln.score >= 10) {
			q10a++;
		}
	}
	printf("Number of confidently mapped reads Q10 %u\n", q10a);
	return;
#endif

	int valid_hash = 0;
	int n_collected = 0;
	int n_proc_contigs = 0;
	int processed_true = 0;
	int bucketed_true = 0;
	int confident = 0;
	int acc = 0;
	int best_hits = 0;
	int true_pos_hits = 0;
	int score = 0;
	int max_votes_inl = 0;
	int max_votes_all = 0;
	int q10 = 0;
	int q30 = 0;
	int q30acc = 0;
	int q10acc = 0;
	int q30processed_true = 0;
	int q30bucketed_true = 0;
	int q10bucketed_true = 0;
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->valid_minhash_f && !r->valid_minhash_rc) continue;
		valid_hash++;
		n_proc_contigs += r->n_proc_contigs;
		true_pos_hits += r->true_n_bucket_hits;
		if(r->collected_true_hit) {
			n_collected++;
		}
		if(r->processed_true_hit) {
			processed_true++;
		}
		if(r->bucketed_true_hit) {
			bucketed_true++;
		}
		if(r->top_aln.score < 5){
			continue;
		}
		confident++;
		if(!r->top_aln.rc) {
			if(pos_in_range(r->ref_pos_l, r->top_aln.ref_start, 20)) {
				r->dp_hit_acc = 1;
			}
		} else {
			if(pos_in_range(r->ref_pos_r, r->top_aln.ref_start, 20)) {
				r->dp_hit_acc = 1;
			}
		}
		acc += r->dp_hit_acc;
		best_hits += r->best_n_bucket_hits;
		score += r->top_aln.score;
		max_votes_inl += r->top_aln.inlier_votes;
		max_votes_all += r->top_aln.total_votes;

		if(r->top_aln.score >= 30) {
			q30++;
			if(r->dp_hit_acc) {
				q30acc++;
			}
			if(r->processed_true_hit) {
				q30processed_true++;
			}
			if(r->bucketed_true_hit) {
				q30bucketed_true++;
			}
		}
		if(r->top_aln.score >= 10) {
			q10++;
			if(r->dp_hit_acc) {
				q10acc++;
			}
			if(r->bucketed_true_hit) {
				q10bucketed_true++;
			}
		}
	}
	float avg_n_proc_contigs = (float) n_proc_contigs/valid_hash;
	float avg_true_pos_hits = (float) true_pos_hits/valid_hash;
	double stddev_true_pos_hits = 0;
	for(uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		if(!r->valid_minhash_f && !r->valid_minhash_rc) continue;
		stddev_true_pos_hits += pow((double)avg_true_pos_hits - r->true_n_bucket_hits, (double) 2);
	}
	stddev_true_pos_hits = sqrt((double)stddev_true_pos_hits/valid_hash);

	printf("Number of reads with valid F or RC hash %u \n", valid_hash);
	printf("Number of mapped reads COLLECTED true hit %u \n", n_collected);
	printf("Number of mapped reads PROC true hit %u \n", processed_true);
	printf("Number of mapped reads BUCK true hit %u \n", bucketed_true);
	printf("Number of confidently mapped reads > 0 %u / accurate %u (%f pct)\n", confident, acc, (float)acc/(float)confident);
	printf("Number of confidently mapped reads Q10 %u / accurate %u (%f pct)\n", q10, q10acc, 100 - 100 * (float)q10acc/(float)q10);
	printf("Number of confidently mapped reads Q30 %u / accurate %u (%f pct)\n", q30, q30acc, (float)q30acc/(float)q30);
	printf("Number of confidently mapped reads Q30 PROCESSED true %u \n", q30processed_true);
	printf("Number of confidently mapped reads Q30 BUCKET true %u \n", q30bucketed_true);
	printf("Number of confidently mapped reads Q10 BUCKET true %u \n", q10bucketed_true);
	printf("Avg number of max bucket hits per read %.8f \n", (float) best_hits/confident);
	printf("MEAN number of bucket hits with TRUE pos per read %.8f \n", (float) avg_true_pos_hits);
	printf("STDDEV number of bucket hits with TRUE pos per read %.8f \n", (float) stddev_true_pos_hits);
	printf("Avg score per read %.8f \n", (float) score/confident);
	printf("Avg max inlier votes per read %.8f \n", (float) max_votes_inl/confident);
	printf("Avg max all votes per read %.8f \n", (float) max_votes_all/confident);
	printf("AVG N_PROC_CONTIGS %.8f \n", avg_n_proc_contigs);

	/*std::string fname("process_contigs");
	fname += std::string("_h");
	fname += std::to_string(params->h);
	fname += std::string("_T");
	fname += std::to_string(params->n_tables);
	fname += std::string("_b");
	fname += std::to_string(params->sketch_proj_len);
	fname += std::string("_w");
	fname += std::to_string(params->ref_window_size);
	fname += std::string("_p");
	fname += std::to_string(params->n_buckets_pow2);
	fname += std::string("_k");
	fname += std::to_string(params->k);
	fname += std::string("_H");
	fname += std::to_string(params->max_count);
	std::ofstream stats_file;
	stats_file.open(fname.c_str(), std::ios::out | std::ios::app);
	if (!stats_file.is_open()) {
		printf("store_ref_idx: Cannot open the CONTIG STATS file %s!\n", fname.c_str());
		exit(1);
	}
	//stats_file << "length,count,m\n";
	for(uint32 i = 0; i < reads.reads.size(); i++) {
				read_t* r = &reads.reads[i];
				if(!r->valid_minhash && !r->valid_minhash_rc) continue;
		stats_file << r->len << "," << r->n_proc_contigs << "," << params->min_n_hits << "\n";
	}
		stats_file.close();
	*/
}