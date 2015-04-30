#include "sam.h"

void print_aln2sam(FILE* samFile, read_t* r, const ref_t& ref);

void store_alns_sam(reads_t& reads, const ref_t& ref, const index_params_t* params) {
	std::string samFname(reads.fname);
	samFname += std::string(".sam");

	FILE* samFile = (FILE*) fopen(samFname.c_str(), "w");
	if (samFile == NULL) {
		printf("alns2sam: Cannot open SAM file: %s!\n", samFname.c_str());
		exit(1);
	}

	for (uint32 i = 0; i < reads.reads.size(); i++) {
		read_t* r = &reads.reads[i];
		print_aln2sam(samFile, r, ref);
	}
	fclose(samFile);
}

#define SAM_FSU   4 // self-unmapped
#define SAM_FSR  16 // self on the reverse strand
void print_aln2sam(FILE* samFile, read_t* r, const ref_t& ref) {
	int flag = 0; // FLAG
	if(r->top_aln.ref_start != 0) {
		seq_t aln_pos = r->top_aln.ref_start;
		if(ref.subsequence_offsets.size() > 1) {
			aln_pos -= ref.subsequence_offsets[r->seq_id];
		}

		r->seq_id = ref.subsequence_offsets.size()-1;
		for(uint32 i = 0; i < ref.subsequence_offsets.size()-1; i++) {
			if(r->top_aln.ref_start >= ref.subsequence_offsets[i] && r->top_aln.ref_start < ref.subsequence_offsets[i+1]) {
				r->seq_id = i;
				break;
			}
		}

		if (r->top_aln.rc) flag |= SAM_FSR;

		// QNAME, FLAG, RNAME
		fprintf(samFile, "%s\t%d\tREF_NAME\t", r->name.c_str(), flag);
		// POS (1-based), MAPQ
		fprintf(samFile, "%d\t%d\t", (int)(aln_pos+1), r->top_aln.score);

		// CIGAR
		fprintf(samFile, "CIGAR");

		// RNEXT, PNEXT, TLEN (print void mate position and coordinate)
		fprintf(samFile, "\t*\t0\t0\t");

		// SEQ, QUAL (print sequence and quality)
		const char* seq = r->top_aln.rc ? r->rc.c_str() : r->seq.c_str();
		for (uint32 i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", "AGCTN"[(int)seq[i]]);
		}
		fprintf(samFile, "\t");
		for (uint32 i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", '*');
		}
		fprintf(samFile, "\n");
	} else { // unmapped read
		int flag = SAM_FSU;
		// QNAME, FLAG
		fprintf(samFile, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", r->name.c_str(), flag);

		// SEQ, QUAL (print sequence and quality)
		for (uint32 i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", "AGCTN"[(int)r->seq[i]]);
		}
		fprintf(samFile, "\t");
		for (uint32 i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", '*');
		}
		fprintf(samFile, "\n");
	}
}
