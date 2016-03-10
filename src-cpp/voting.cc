#include <algorithm>
#include "voting.h"

inline bool is_unique_kmer(const std::vector<kmer_enc_t>& sorted_kmers, const uint32 idx) {
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

void record_vote_discretized(const kmer_enc_t& rkmer, const kmer_enc_t& ckmer, const int max_rpos, const int max_cpos, const int offset, std::vector<int>& votes) {
	const int bin_size = params->k2;
	
	const int rbin_idx = rkmer.pos;
	const int rrange_start = rbin_idx*bin_size;
	int rrange_end = rrange_start + bin_size;	
	if(rrange_end > max_rpos) rrange_end = max_rpos;
	
	const int cbin_idx = ckmer.pos;			
	const int crange_start = cbin_idx*bin_size;
	int crange_end = crange_start + bin_size;
	if(crange_end > max_cpos) crange_end = max_cpos;	

	const int s = crange_start - rrange_end + 1;
	const int t = crange_end - rrange_start;

	for(int k = s; k < t; k++) {
		votes[offset + k]++;
	}
}

// finds the matching contig and read kmers and records their votes
// returns true if there was at least one match between the kmers
bool vote_cast_and_count(const std::vector<kmer_enc_t>& read_kmers, const std::vector<kmer_enc_t>& contig_kmers, const int max_rpos, const int max_cpos, std::vector<int>& votes) {
	votes.resize(max_rpos + max_cpos); // minimal req: last read kmer + first contig kmer
	const int cstart_pos = max_rpos;
	uint32 skip = 0;
	bool any_matches = false;
	for(int i = 0; i < read_kmers.size(); i++) {
		if(read_kmers[i].hash == 0) continue;
		if(!is_unique_kmer(read_kmers, i)) continue;
		for(uint32 j = skip; j < contig_kmers.size(); j++) {
			if(!is_unique_kmer(contig_kmers, j)) continue;
			if(read_kmers[i].hash == contig_kmers[j].hash) {
				any_matches = true;
				record_vote_discretized(read_kmers[i], contig_kmers[j], max_rpos, max_cpos,  cstart_pos, votes);
			} else if(contig_kmers[j].hash > read_kmers[i].hash) {
					skip = j;
					break;
			}
		}
	}
	return any_matches;
}

void find_kmer_matches(const std::vector<kmer_enc_t>& read_kmers, const std::vector<kmer_enc_t>& contig_kmers, std::vector<kmer_match_t>& matches) {
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
void prefsum(const std::vector<int>& votes, std::vector<int> votes_prefsum) {
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
void find_max_vote_mid(const std::vector<int> votes_prefsum, int& max_val, int& max_idx, const int conv_range, const bool selective, const int exclude_pos, const int exclude_range) {
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
void prepare_voting_read_kmers(const kmer_cipher_t* rkmers, std::vector<kmer_enc_t>& rvk) {
	const int bin_size = params->k2;
	for(int i = 0; i < rvk.size(); i++) {
		rvk[i] = kmer_enc_t(rkmers[i], i/bin_size);
	}
	std::sort(rvk.begin(), rvk.end(), kmer_enc_t::comp());
}

// sparse, pos scrambled within bin_size
void prepare_voting_contig_kmers(const kmer_cipher_t* ckmers, std::vector<kmer_enc_t>& cvk) {
	const int bin_size = params->k2;
	const int bin_size_sampled = bin_size/params->sampling_intv;
	for(int i = 0; i < cvk.size(); i++) {
		cvk[i] = kmer_enc_t(ckmers[i], i/bin_size_sampled);
	}
	std::sort(cvk.begin(), cvk.end(), kmer_enc_t::comp());
}

void voting_task::process(voting_results& out) {
	const int n_read_kmers = get_read_data_len();
	std::vector<kmer_enc_t> rvk(n_read_kmers);
	prepare_voting_read_kmers(get_read(), rvk);
	int n_contigs = get_n_contigs();
	for(int i = 0; i < n_contigs; i++) {
		const int n_contig_kmers = get_contig_data_len(i);
		std::vector<kmer_enc_t> cvk(n_contig_kmers);
		prepare_voting_contig_kmers(get_contig(i), cvk);
		
		// find matches and count votes
		std::vector<int> votes;
		bool any_votes = vote_cast_and_count(rvk, cvk, n_read_kmers, get_n_kmers(get_contig_len(i), params->k2), votes);
		if(any_votes) return;
	
		// convolution
		std::vector<int> votes_prefsum;
		prefsum(votes, votes_prefsum);
		// max vote
		voting_results cur;
		find_max_vote_mid(votes_prefsum , cur.best_score[voting_results::topid::BEST], cur.local_pos[voting_results::topid::BEST], params->delta_inlier, false, 0, 0);
		cur.local_pos[0] -= n_read_kmers;
		// second best
		find_max_vote_mid(votes_prefsum , cur.best_score[voting_results::topid::SECOND], cur.local_pos[voting_results::topid::SECOND], params->delta_inlier, true, cur.local_pos[voting_results::topid::BEST], params->delta_x*params->delta_inlier);
		cur.local_pos[1] -= n_read_kmers;
		out.compare_and_update(cur);
	}
}

// input:
// - tasks of encrypted read + contig kmers
// ouput: 
// - best candidate per task
// - stats
void run_voting(const std::vector<voting_task*>& tasks, std::vector<voting_results>& results, voting_stats& stats) {
	int sum = 0;
	int n_nonzero = 0;
	for(size_t i = 0; i < tasks.size(); i++) {
		tasks[i]->process(results[i]);
		if(results[i].best_score[voting_results::topid::BEST] > 0) {
			sum += results[i].best_score[voting_results::topid::BEST];
			n_nonzero++;
		}
	}
	if(n_nonzero > 0) {
		stats.avg_score = sum/n_nonzero;
	}
}