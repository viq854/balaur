#include "voting.h"

static const int N_INIT_ANCHORS_MAX = (getenv("N_INIT_ANCHORS_MAX") ? atoi(getenv("N_INIT_ANCHORS_MAX")) : 20);
#define UNIQUE1 0
#define UNIQUE2 1
void ransac(std::vector<kmer_enc_t>& read_ciphers, std::vector<kmer_enc_t>& contig_ciphers, int* n_votes, int* pos) {
	int anchors[2][N_INIT_ANCHORS_MAX];
	uint32 n_anchors[2] = { 0 };
	int anchor_median_pos[2] = { 0 };

	for(int p = 0; p < 3; p++) {
		// start from the beginning in each pass
		uint32 idx_q = 0;
		uint32 idx_r = 0;
		while(idx_q < read_ciphers.size() && idx_r < contig_ciphers.size()) {
			if(contig_ciphers[idx_r].hash == read_ciphers[idx_q].hash) { // MATCH
				int match_aln_pos = contig_ciphers[idx_r].pos - read_ciphers[idx_q].pos;
				if(p == 0) { // -------- FIRST PASS: collect unique matches ---------
					if(is_unique_kmer(read_ciphers, idx_q)) {
						if(n_anchors[UNIQUE1] < params->n_init_anchors) {
							anchors[UNIQUE1][n_anchors[UNIQUE1]] = match_aln_pos;
							n_anchors[UNIQUE1]++;
						} else { // found all the anchors
							break;
						}
					}
					idx_q++;
					idx_r++;
			} else if(p == 1) { // --------- SECOND PASS: collect DIFF unique matches -------
					if(n_anchors[UNIQUE1] == 0) break; // no unique matches exist
					// compute the median first pass position
					if(anchor_median_pos[UNIQUE1] == 0) {
						std::sort(&anchors[UNIQUE1][0], &anchors[UNIQUE1][0] + n_anchors[UNIQUE1]);
						anchor_median_pos[UNIQUE1] = anchors[UNIQUE1][(n_anchors[UNIQUE1]-1)/2];
					}
					if(is_unique_kmer(read_ciphers, idx_q)) {
						// if this position is not close to the first pick
						if(!pos_in_range_sig(match_aln_pos, anchor_median_pos[UNIQUE1], params->delta_x*params->delta_inlier)) {
							if(n_anchors[UNIQUE2] < params->n_init_anchors) {
								anchors[UNIQUE2][n_anchors[UNIQUE2]] = match_aln_pos;
								n_anchors[UNIQUE2]++;
							} else { // found all the anchors
								break;
							}
						}
					}
					idx_q++;
					idx_r++;
				} else { // --------- THIRD PASS: count inliers to each candidate position -------
					if(n_anchors[UNIQUE1] == 0) break;
					// find the median alignment positions
					if(n_anchors[UNIQUE2] > 0 && anchor_median_pos[UNIQUE2] == 0) {
						std::sort(&anchors[UNIQUE2][0], &anchors[UNIQUE2][0] + n_anchors[UNIQUE2]);
						anchor_median_pos[UNIQUE2] = anchors[UNIQUE2][(n_anchors[UNIQUE2]-1)/2];
					}
					// count inliers
					for(int i = 0; i < 2; i++) {
						if(n_anchors[i] == 0) continue;
						if(pos_in_range_sig(anchor_median_pos[i], match_aln_pos, params->delta_inlier)) {
							n_votes[i]++;
							pos[i] += match_aln_pos;
						}
					}
					if(idx_r < (contig_ciphers.size()-1) && contig_ciphers[idx_r + 1].hash == contig_ciphers[idx_r].hash) {
						// if the next kmer is still a match in the reference, check it in next iteration as an inlier
						idx_r++;
					} else {
						idx_q++;
						idx_r++;
					}
				}
			} else if(read_ciphers[idx_q].hash < contig_ciphers[idx_r].hash) {
				idx_q++;
			} else {
				idx_r++;
			}
		}
	}
	for(int i = 0; i < 2; i++) {
		if(n_votes[i] > 0) {
			pos[i] = pos[i]/(int)n_votes[i];
		}
	}

	if(n_anchors[UNIQUE1] < 5) {
		n_votes[0] = 0;
		n_votes[1] = 0;
	}
}