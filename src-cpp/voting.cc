#include <algorithm>
#include "voting.h"

inline bool is_unique_kmer(const VectorKmerEnc& sorted_kmers, const uint32 idx) {
	bool unique = true;
	kmer_cipher_t hash = sorted_kmers[idx].hash;
	if(idx > 0 && sorted_kmers[idx-1].hash == hash) {
		unique = false;
	} 
	if(idx < sorted_kmers.size()-1 && sorted_kmers[idx+1].hash == hash) {
		unique = false;
	}

	return unique;
}

inline void record_vote(const kmer_enc_t& rkmer, const kmer_enc_t& ckmer, const int offset, std::vector<int>& votes) {
	int match_aln_pos = ckmer.pos - rkmer.pos;
	votes[offset + match_aln_pos]++;
}


void get_vote_pos_range(const kmer_enc_t& rkmer, const kmer_enc_t& ckmer, const int max_rpos, const int max_cpos, std::pair<int, int>& range) {
	const int bin_size = params->bin_size;
	
	const int rbin_idx = rkmer.pos;
	const int rrange_start = rbin_idx*bin_size;
	int rrange_end = rrange_start + bin_size;	
	if(rrange_end > max_rpos) rrange_end = max_rpos;
	
	const int cbin_idx = ckmer.pos;			
	const int crange_start = cbin_idx*bin_size;
	int crange_end = crange_start + bin_size;
	if(crange_end > max_cpos) crange_end = max_cpos;	

	const int s = crange_start - rrange_end; // + 1;
	const int t = crange_end - rrange_start;
	range.first = s;
	range.second = t;
}

void record_vote_discretized(const kmer_enc_t& rkmer, const kmer_enc_t& ckmer, const int max_rpos, const int max_cpos, const int offset, std::vector<int>& votes) {
	std::pair<int, int> bin_range;
	get_vote_pos_range(rkmer, ckmer, max_rpos, max_cpos, bin_range);
	const int s = bin_range.first;
	const int t = bin_range.second;
	for(int k = s; k < t; k++) {
		votes[offset + k]++;
	}
}

// finds the matching contig and read kmers and records their votes
// returns true if there was at least one match between the kmers
bool vote_cast_and_count(const VectorKmerEnc& read_kmers, const VectorKmerEnc& contig_kmers, const int max_rpos, const int max_cpos, std::vector<int>& votes) {
	votes.resize(max_rpos + max_cpos); // minimal req: last read kmer + first contig kmer
	const int cstart_pos = max_rpos;
	uint32 skip = 0;
	bool any_matches = false;
	for(int i = 0; i < read_kmers.size(); i++) {
		if(read_kmers[i].hash == 0) continue;
		if(!is_unique_kmer(read_kmers, i)) continue;
		//uint32 j;
		for(uint32 j = skip; j < contig_kmers.size(); j++) {
			if(!is_unique_kmer(contig_kmers, j)) continue;
			if(read_kmers[i].hash == contig_kmers[j].hash) {
				any_matches = true;
				record_vote_discretized(read_kmers[i], contig_kmers[j], max_rpos, max_cpos, cstart_pos, votes);
			} else if(contig_kmers[j].hash > read_kmers[i].hash) {
					skip = j;
					break;
			}
		}
		//if(j == contig_kmers.size()) break; 
	}
	return any_matches;
}

void collect_kmer_matched_for_pos(const int aln_pos, const VectorKmerEnc& read_kmers, const VectorKmerEnc& contig_kmers, 
const int max_rpos, const int max_cpos, std::vector<cipher_match_t>& matches) {

	uint32 skip = 0;
        for(int i = 0; i < read_kmers.size(); i++) {
                if(read_kmers[i].hash == 0) continue;
                if(!is_unique_kmer(read_kmers, i)) continue;
                for(uint32 j = skip; j < contig_kmers.size(); j++) {
                        if(!is_unique_kmer(contig_kmers, j)) continue;
                        if(read_kmers[i].hash == contig_kmers[j].hash) {
				std::pair<int, int> bin_range;
				get_vote_pos_range(read_kmers[i], contig_kmers[j], max_rpos, max_cpos, bin_range);
				if(aln_pos >= bin_range.first && aln_pos <= bin_range.second) {
					cipher_match_t m;
					m.rk = read_kmers[i];
					m.ck = contig_kmers[j];
					matches.push_back(m);
				}
                        } else if(contig_kmers[j].hash > read_kmers[i].hash) {
                                        skip = j;
                                        break;
                        }
                }
        }
	//std::sort(matches.begin(), matches.end(), cipher_match_t::comp());

	//if(matches.size() == 0) {
        //        std::cout << "ERROR: could not find a match for pos\n";
        //        exit(-1);
        //}

	// testing...	
	if(matches.size() > 0) {
		int min = 10000000;
        	int min_idx = 0;
        	for(int i = 0; i < matches.size(); i++) {
                	if(matches[i].rk.orig_pos < min) {
                        	min = matches[i].rk.orig_pos;
                        	min_idx = i;
                	}
        	}
		cipher_match_t tmp = matches[0];
		matches[0] = matches[min_idx];
		matches[min_idx] = tmp;
	}
}


void find_kmer_matches(const VectorKmerEnc& read_kmers, const VectorKmerEnc& contig_kmers, std::vector<kmer_match_t>& matches) {
	uint32 skip = 0;
	for(uint32 i = 0; i < read_kmers.size(); i++) {
		if(read_kmers[i].hash == 0) continue;
		if(!is_unique_kmer(read_kmers, i)) continue;
		for(uint32 j = skip; j < contig_kmers.size(); j++) {
			if(!is_unique_kmer(contig_kmers, j)) continue;
			if(read_kmers[i].hash == contig_kmers[j].hash) {
				matches.push_back(kmer_match_t(contig_kmers[j].pos, read_kmers[i].pos));
			} else if(contig_kmers[j].hash > read_kmers[i].hash) {
				skip = j;
				break;
			}
		}
	}
}

/*void vote_cast_and_count_chaining(const seq_t rlen, const seq_t clen,  const std::vector<kmer_enc_t>& read_kmers, const std::vector<kmer_enc_t>& contig_kmers, std::vector<int>& votes) {
	std::vector<kmer_match_t> matches;
	find_kmer_matches(read_kmers, contig_kmers, matches);
	if(matches.size() == 0) return;
	std::sort(matches.begin(), matches.end());

	// find maximal overlapping fragments
	std::vector<int> fragment_weights;
	std::vector<kmer_match_t> fragments;
	int idx = 0;
	fragments.push_back(kmer_match_t(matches[0].cpos, matches[0].rpos));
	fragment_weights.push_back(params->k2);
	for(int i = 1; i < matches.size(); i++) {
		int match_aln_pos = matches[i].cpos - matches[i].rpos;
		int delta_contig = matches[i].cpos - matches[i-1].cpos;		
		int delta_read = matches[i].rpos - matches[i-1].rpos;
		int last_pos = fragments[idx].cpos + fragment_weights[idx];
		if(delta_contig == delta_read && matches[i].cpos < (last_pos + 1)) {
			fragment_weights[idx] += matches[i].cpos + params->k2 - last_pos;
		} else {
			int r = 0;
			if(matches[i].rpos > matches[i-1].rpos) { // indels
				int d1 = last_pos - matches[i].cpos;
				int d2 = fragments[idx].rpos + fragment_weights[idx] - matches[i].rpos;
				if(d1 > 0) r = d1;
				if(d2 > 0) r = d1 > d2 ? d1 : d2;	
			}
			fragments.push_back(kmer_match_t(matches[i].cpos + r, matches[i].rpos + r));
			fragment_weights.push_back(params->k2 - r);
			idx++;
		}
	}
	std::vector<int> best_chain;
	chain_fragments_dp(fragments, fragment_weights, best_chain);

	votes.resize(clen + rlen);
	const int cstart_pos = rlen;
	for(int i = 0; i < best_chain.size(); i++) {      
        	int match_aln_pos = fragments[best_chain[i]].cpos - fragments[best_chain[i]].rpos;
		votes[cstart_pos + match_aln_pos] += get_n_kmers(fragment_weights[best_chain[i]], params->k2);
	}
	//std::vector<int> votes(ref_contig.len + rlen);
	//for(int i = 1; i < seqan::length(global_chain)-1; i++) {	
	//	int match_aln_pos = seqan::leftPosition(global_chain[i], 0) - seqan::leftPosition(global_chain[i], 1);
	//	votes[cstart_pos + match_aln_pos] += 20*(seqan::rightPosition(global_chain[i], 0) - seqan::leftPosition(global_chain[i], 0));
	//}
}*/

// optimized implementation of the convolution step using the prefix sum
void prefsum(const std::vector<int>& votes, std::vector<int>& votes_prefsum) {
	votes_prefsum.resize(votes.size());
	votes_prefsum[0] = votes[0];
	for(int i = 1; i < votes.size(); i++) {
		votes_prefsum[i] = votes[i] + votes_prefsum[i-1];
	}
}

// find the max votes value and position
// in the convolved votes array (optimization: convolution is done on the fly using the prefix sum)
// returns the middle position in the range of positions with the same max value
// can selectively exclude the vecitiny of a given position from consideration (e.g. used for finding second best)
void find_max_vote_mid(const std::vector<int>& votes_prefsum, int& max_val, int& max_idx, const int conv_range, const bool selective, const int exclude_pos, const int exclude_range) {
	max_val = 0;
	for(int i = 0; i < votes_prefsum.size(); i++) {
		if(selective && pos_in_range(exclude_pos, i, exclude_range)) continue;
		uint32 start = i > conv_range ? i - conv_range - 1 : 0;
		uint32 end = i + conv_range >=  votes_prefsum.size() ? votes_prefsum.size() - 1: i + conv_range;
		int votes_conv = votes_prefsum[end] - votes_prefsum[start];
		if(votes_conv > max_val) {
			max_val = votes_conv;
			max_idx = i;
		}
	}
	if(selective) return; // && max_val == 0) return;
	// pick the middle position in the max-val range
	int i = max_idx; 
	while(i < votes_prefsum.size()) {
		uint32 start = i > conv_range ? i - conv_range - 1 : 0;
		uint32 end = i + conv_range >=  votes_prefsum.size() ? votes_prefsum.size() - 1: i + conv_range;
		int votes_conv = votes_prefsum[end] - votes_prefsum[start];
		if(votes_conv == max_val) {
			i++;
		} else {
			break;
		}
	}
	max_idx = (i + max_idx)/2;
}

// dense, pos scrambled within bin_size
void prepare_voting_read_kmers(const kmer_cipher_t* rkmers, VectorKmerEnc& rvk) {
	for(int i = 0; i < rvk.size(); i++) {
		rvk[i] = kmer_enc_t(rkmers[i], i/params->bin_size);
		rvk[i].orig_pos = i; // temp
	}
	std::sort(rvk.begin(), rvk.end(), kmer_enc_t::comp());
}

// sparse, pos scrambled within bin_size
void prepare_voting_contig_kmers(const kmer_cipher_t* ckmers, VectorKmerEnc& cvk) {
	int bin_size_sampled = params->bin_sampling() ? params->bin_size/params->sampling_intv : 1;
	for(int i = 0; i < cvk.size(); i++) {
		int pos = params->bin_sampling() ? i/bin_size_sampled : i*params->sampling_intv;
		cvk[i] = kmer_enc_t(ckmers[i], pos);
	}
	std::sort(cvk.begin(), cvk.end(), kmer_enc_t::comp());
}

void voting_task::process(voting_results& out) {
	const int n_read_kmers = get_read_data_len();
	//std::vector<kmer_enc_t> rvk(n_read_kmers);
	if(n_read_kmers < 0) {
		std::cout << "ERROR: n_read_kmers < 0\n";
		exit(-1); 
	}
	out.rvk.resize(n_read_kmers);
	prepare_voting_read_kmers(get_read(), out.rvk);
	int n_contigs = get_n_contigs();

	for(int i = 0; i < n_contigs; i++) {
		voting_results cur;
		const int n_contig_kmers = get_contig_data_len(i);
		//std::vector<kmer_enc_t> cvk(n_contig_kmers);
		cur.cvk.resize(n_contig_kmers);
		prepare_voting_contig_kmers(get_contig(i), cur.cvk);
		cur.contig_len = get_contig_len(i);		

		// find matches and count votes
		std::vector<int> votes;
		bool any_votes = vote_cast_and_count(out.rvk, cur.cvk, n_read_kmers, get_n_kmers(get_contig_len(i), params->k2), votes);
		if(!any_votes) {
			continue;
		}
		
		// convolution
		std::vector<int> votes_prefsum;
		prefsum(votes, votes_prefsum);

		// max vote
		//voting_results cur;
		find_max_vote_mid(votes_prefsum, cur.best_score[voting_results::topid::BEST], cur.local_pos[voting_results::topid::BEST], params->delta_inlier, false, 0, 0);

		// second best
		find_max_vote_mid(votes_prefsum, cur.best_score[voting_results::topid::SECOND], cur.local_pos[voting_results::topid::SECOND], params->delta_inlier, true, cur.local_pos[voting_results::topid::BEST], 50); //params->delta_x);
		cur.local_pos[0] -= n_read_kmers;
		cur.local_pos[1] -= n_read_kmers;
		cur.contig_id[0] = contig_ids[i];
		cur.contig_id[1] = contig_ids[i];
		cur.global_pos[0] = global_pos[i];
		cur.global_pos[1] = global_pos[i];

		out.compare_and_update(cur);
	}
	//if(out.best_score[voting_results::topid::BEST] > 0) std::cout << out.rid << " best score " << out.best_score[voting_results::topid::BEST] << " " << out.best_score[voting_results::topid::SECOND] << "\n";
}

// input:
// - tasks of encrypted read + contig kmers
// ouput: 
// - best candidate per task
// - stats
void run_voting(const std::vector<voting_task*>& tasks, std::vector<voting_results>& results, voting_stats& stats) {
	//int sum = 0;
	//int n_nonzero = 0;
	
	omp_set_num_threads(params->n_threads);
	std::cout << "Running the voting task on " << params->n_threads << " threads\n";
	double start_time = omp_get_wtime(); 
	#pragma omp parallel for //reduction(+:sum, n_nonzero)
	for(size_t i = 0; i < tasks.size(); i++) {
		results[i].rid = tasks[i]->rid;
		results[i].rc = tasks[i]->strand;
		results[i].key1_xor_pad = tasks[i]->key1_xor_pad; //temp
        	results[i].key2_mult_pad = tasks[i]->key2_mult_pad;
		tasks[i]->process(results[i]);
		if(results[i].best_score[voting_results::topid::BEST] > 0) {
			//sum += results[i].best_score[voting_results::topid::BEST];
			//n_nonzero++;
			// collect matching kmers
			if(results[i].rvk.size() == 0 || results[i].cvk.size() == 0) {
				std::cout << "ERROR: empty cipher vectors \n";
				exit(-1); 
			}
			const int aln_pos = results[i].local_pos[voting_results::topid::BEST];
			const int cid = results[i].contig_id[voting_results::topid::BEST];
			collect_kmer_matched_for_pos(aln_pos, results[i].rvk, results[i].cvk, tasks[i]->get_read_data_len(), get_n_kmers(results[i].contig_len, params->k2),
				results[i].kmer_matches); 
			//if(results[i].kmer_matches.size() == 0) {
                        //        std::cout << "ERROR: did not collect any matches \n";
                        //        std::cout << aln_pos << " " << cid << "\n";
			//	exit(1);
			//}
				/*uint32 skip = 0;	
				for(int x = 0; x < results[i].rvk.size(); x++) {
					if(!is_unique_kmer(results[i].rvk, x)) continue;
					for(uint32 j = skip; j < results[i].cvk.size(); j++) {
						if(!is_unique_kmer(results[i].cvk, j)) continue;
						if(results[i].rvk[x].hash == results[i].cvk[j].hash) {
							std::cout << "FOUND MATCH " << results[i].rvk[x].pos << " " << results[i].cvk[j].pos << "; ";
							std::pair<int, int> bin_range;
							get_vote_pos_range(results[i].rvk[x], results[i].cvk[j], tasks[i]->get_read_data_len(), get_n_kmers(results[i].contig_len, params->k2), 
								bin_range);
							
							int max_rpos = tasks[i]->get_read_data_len();
							int max_cpos = get_n_kmers(results[i].contig_len, params->k2);
							const int rbin_idx = results[i].rvk[x].pos;
        						const int rrange_start = rbin_idx*params->bin_size;
        						int rrange_end = rrange_start + params->bin_size;
        						if(rrange_end > max_rpos) rrange_end = max_rpos;

        						const int cbin_idx = results[i].cvk[j].pos;
        						const int crange_start = cbin_idx*params->bin_size;
        						int crange_end = crange_start + params->bin_size;
        						if(crange_end > max_cpos) crange_end = max_cpos;

        						const int s = crange_start - rrange_end; // + 1;
        						const int t = crange_end - rrange_start;
							std::cout << " RANGE: " << bin_range.first << " " << bin_range.second << "; " << s << " " << t << " " << max_rpos << " " << max_cpos << "\n";
						} else if(results[i].cvk[j].hash > results[i].rvk[x].hash) {
							skip = j;
							break;
						}
					}
				}
				exit(-1);
				*/
                        //}
		}
	}
	//if(n_nonzero > 0) {
	//	stats.avg_score = sum/n_nonzero;
	//}
	printf("Voting runtime : %.2f sec\n", omp_get_wtime() - start_time);
}
