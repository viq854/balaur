#if(SEQAN_LIB)
#include <seqan/seeds.h>
#include <seqan/graph_algorithms.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#endif

#if(CLASP_LIB)
extern "C" {
#include "clasp_v1_1/libs/sltypes.h"
#include "clasp_v1_1/libs/slchain.h"
#include "clasp_v1_1/libs/container.h"
#include "clasp_v1_1/libs/debug.h"
#include "clasp_v1_1/libs/mathematics.h"
#include "clasp_v1_1/libs/sort.h"
#include "clasp_v1_1/libs/vqueue.h"
#include "clasp_v1_1/libs/vebtree.h"
#include "clasp_v1_1/libs/bintree.h"
#include "clasp_v1_1/libs/rangetree.h"
}
#endif

#include "lemon/list_graph.h"
#include <lemon/connectivity.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/graph_components.h>

#define CHAIN_VOTES 1
using namespace lemon;
struct node_comp {
	ListDigraph::NodeMap<int>* order;
	node_comp(ListDigraph::NodeMap<int>* _order) {
		 order = _order;
	}
	bool operator()( const ListDigraph::Node& x, const ListDigraph::Node& y ) const {
		return (*order)[x] < (*order)[y];
	}
};

void chain_fragments_dp(std::vector<kmer_match_t>& fragments, std::vector<int>& fragment_weights, std::vector<int>& best_chain) {
	int n_frags = fragments.size();	
	
	// initialize DAG
	ListDigraph g;
	ListDigraph::NodeMap<int> frag_weights(g, 0); 
	ListDigraph::ArcMap<double> gap_weights(g, 0);
	ListDigraph::NodeMap<int> frag_ids(g);
	std::map<int, ListDigraph::Node> fragId_to_node;

	// add source node
	ListDigraph::Node source = g.addNode();
	ListDigraph::Node sink = g.addNode();
	frag_weights.set(source, 0);
	frag_weights.set(sink, 0);

	// each fragment is a node
	for(int i = 0; i < n_frags; i++) {
		// create new node with weight
		ListDigraph::Node u = g.addNode();
		frag_weights.set(u, 5*fragment_weights[i]);
		// link to source and sink with gap penalty 0
		ListDigraph::Arc a = g.addArc(source, u);
		ListDigraph::Arc b = g.addArc(u, sink);
		gap_weights.set(a, 0);
		gap_weights.set(b, 0);
		frag_ids.set(u, i);
		fragId_to_node.insert(std::make_pair(i, u));
	}
	int n_edges = 0;
	int mm_score = 0;
	int gap_score = 10;
	int fixed_cost = 2;
	for(int i = 0; i < n_frags; i++) {
		ListDigraph::Node u = fragId_to_node.find(i)->second;
		for(int j = i+1; j < n_frags; j++) {
			int delta_x = fragments[j].cpos - (fragments[i].cpos + fragment_weights[i]);
			int delta_y = fragments[j].rpos - (fragments[i].rpos + fragment_weights[i]);
			if(delta_x >= 0 && delta_y >= 0) {
				// add edge
				ListDigraph::Node v = fragId_to_node.find(j)->second;
				ListDigraph::Arc  a = g.addArc(u, v);
				
				int gap = delta_x > delta_y ? delta_x - delta_y : delta_y - delta_x;
				int mm = delta_x > delta_y ? delta_y : delta_x;
				gap_weights.set(a, gap*gap_score + mm_score*mm +fixed_cost); // custom gap scoring function
				n_edges++;
			}
		}
	}

	// topological sort
	ListDigraph::NodeMap<int> order(g);
	topologicalSort(g, order);

	std::vector<ListDigraph::Node> nodes;
	for(ListDigraph::NodeIt u(g); u != INVALID; ++u) {
		nodes.push_back(u);
	}
	std::sort(nodes.begin(), nodes.end(), node_comp(&order));
	ListDigraph::NodeMap<int> best_score(g, 0);
	best_score.set(source, 0);
	ListDigraph::NodeMap<ListDigraph::Node> pred(g);
	for(int i = 1; i < n_frags + 2; i++) {
		ListDigraph::Node u = nodes[i]; // ordered	
		
		int max_score = 0;
		ListDigraph::Node pred_n;
		// get incoming edges
		for(ListDigraph::InArcIt e(g,u); e != INVALID; ++e ) {
			ListDigraph::Node v = g.source(e);
			// get the best score to this node
			int score_v = best_score[v];
			int score = score_v - gap_weights[e] + frag_weights[u];
			if(score >= max_score) {
				max_score = score;
				pred_n = v;
			}
		}
		best_score.set(u, max_score);
		pred.set(u, pred_n);
	}
	// start at sink
	ListDigraph::Node n = sink;
	int best_chain_score = best_score[sink];
	while(n != source) {
		n = pred[n];
		if(n != source) best_chain.push_back(frag_ids[n]);
		//std::cout << "frag " << frag_ids[n] << " " << best_score[n] << "\n";
	}
	//std::cout << "best chain score " << best_chain_score << "\n";
}

void chain_fragments_local(std::vector<kmer_match_t>& fragments, std::vector<int>& fragment_weights) {
#if(CLASP_LIB)
	int n_frags = fragments.size();
	slmatch_t* frags =  (slmatch_t *) malloc(n_frags*sizeof(slmatch_t));
	for (int i = 0; i < n_frags; i++) {
		slmatch_t* frag = &frags[i];
		bl_slmatchInit(frag, 0);
		frag->p = fragments[i].cpos;
		frag->q = fragment_weights[i];
		frag->i = fragments[i].rpos;
		frag->j = fragment_weights[i];
		frag->scr = fragment_weights[i]*10;
		frag->idx = i;
		frag->subject = 0;
	}
	int eps = 0;
	int lambda = 1;
	int maxgap = -1;
	bl_slChainSop(frags, n_frags, eps, lambda, maxgap);

	for(int i = 0; i < n_frags; i++){
		slmatch_t *match = &frags[i];
		std::cout << "FRAG " << match->i << " " <<  match->i + match->j - 1 << " " << match->p << " " <<  match->p + match->q - 1 << " " <<  match->scr << " \n";
		if (match->chain) {
			slchain_t *chain = (slchain_t *) match->chain;
			if (chain->scr >= 0 && bl_containerSize(chain->matches) >= 1) {
				fprintf(stderr, "CHAIN %d\t%d\t%d\t%d\t%.3f\n", chain->i, chain->i + chain->j - 1, chain->p, chain->p + chain->q - 1, chain->scr);
				for (int k = 0; k < bl_containerSize(chain->matches); k++){
						slmatch_t *frag = *(slmatch_t **) bl_containerGet(chain->matches, k);
						fprintf(stderr, "%d\t%d\t%d\t%d\t%.3f\n", frag->i, frag->i + frag->j - 1, frag->p, frag->p + frag->q - 1, frag->scr);
				}
			}
			//bl_slchainDestruct(chain);
			//free(chain);
			//match->chain = NULL;
		}
	}
#endif
}

void chain_fragments_global(std::vector<kmer_match_t>& fragments, std::vector<int>& fragment_weights) {
#if(SEQAN_LIB)
	seqan::String<seqan::Seed<int, seqan::MultiSeed>> s_fragments;
        for (int i = 0; i < fragments.size(); i++) {
                seqan::Seed<int, seqan::MultiSeed> seed(2);
                seqan::setLeftPosition(seed, 0, fragments[i].cpos);
                seqan::setRightPosition(seed, 0, fragments[i].cpos + fragment_weights[i]);
                seqan::setLeftPosition(seed, 1, fragments[i].rpos);
                seqan::setRightPosition(seed, 1, fragments[i].rpos + fragment_weights[i]);
                seqan::setWeight(seed, fragment_weights[i]*10);
                seqan::appendValue(s_fragments, seed);
        }

        int gap_diff_score = 3;
        int gap_m_score = 1;
        seqan::String<seqan::Seed<int, seqan::MultiSeed> > global_chain;
        seqan::Score<int, seqan::ChainSoP> gap_scoring(0, gap_m_score, gap_diff_score);
        int chain_score = seqan::globalChaining(s_fragments, global_chain, gap_scoring);

        int p1_x = seqan::leftPosition(global_chain[1], 0);
        int p1_y = seqan::leftPosition(global_chain[1], 1);
        int min_p_start = (p1_x < p1_y) ? p1_x : p1_y;
        int diff_p_start = (p1_x < p1_y) ? p1_y - p1_x : p1_x - p1_y;
        chain_score += min_p_start*gap_m_score;
        chain_score += diff_p_start*gap_diff_score;

        int p2_x = seqan::leftPosition(global_chain[seqan::length(global_chain) - 1], 0) - seqan::leftPosition(global_chain[seqan::length(global_chain) - 2], 0);
        int p2_y = seqan::leftPosition(global_chain[seqan::length(global_chain) - 1], 1) - seqan::leftPosition(global_chain[seqan::length(global_chain) - 2], 1);
        int min_p_end = (p2_x < p2_y) ? p2_x : p2_y;
        int diff_p_end = (p2_x < p2_y) ? p2_y - p2_x : p2_x - p2_y;
        chain_score += min_p_end*gap_m_score;
        chain_score += diff_p_end*gap_diff_score;


        //for (int i = 0; i < seqan::length(global_chain); i++) {
        //      std::cout <<  seqan::leftPosition(global_chain[i], 0) << " " <<  seqan::leftPosition(global_chain[i], 1) << ";";
        //}
        //std::cout << "(CHAIN SCORE =" << chain_score << ")\n";
        //std::cout << chain_score << "\n";

        // HIS over the fragments
        //seqan::String<unsigned int> weights;
        //seqan::resize(weights, fragment_weights.size(), 0);
        //for(int i = 0; i < fragment_weights.size(); i++) {
        //      seqan::assignProperty(weights, i, fragment_weights[i]);
        //}
        //typedef seqan::Position<seqan::String<unsigned int> >::Type TPosition;
        //seqan::String<TPosition> spos;
        //w = heaviestIncreasingSubsequence(fragments, weights, spos);

        //for (int i = 0; i < fragments.size(); i++) {
        //      std::cout << fragments[i] << "(Weight=" << seqan::getProperty(weights, i) << "),";
        //}
        //std::cout << "\n" << "His: \n";
        //for (int i = length(spos) - 1; i >= 0; i--) {
        //      std::cout << fragments[spos[i]] <<  ',';
        //}
        //std::cout << "(HIS Weight=" << w << ")\n";
#endif
}