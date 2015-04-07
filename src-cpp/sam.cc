#include "sam.h"

void print_aln2sam(FILE* samFile, read_t* r, const ref_t& ref);

void store_alns_sam(reads_t& reads, const ref_t& ref, const index_params_t* params) {
	std::string samFname(reads.fname);
	samFname += std::string(".sam");

	FILE* samFile = (FILE*) fopen(samFname.c_str(), "w");
	if (samFile == NULL) {
		printf("alns2sam: Cannot open SAM file: %s!\n", samFname);
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
void print_aln2sam(FILE* samFile, const ref_t& ref, const read_t* r) {
	int flag = 0; // FLAG
	if(r->aln.ref_start != 0) {
		if(ref.subsequence_offsets.size() > 1) {
			r->aln.read_start -= ref.subsequence_offsets[r->seq_id];
		}

		if (r->aln.rc) flag |= SAM_FSR;

		// QNAME, FLAG, RNAME
		fprintf(samFile, "%s\t%d\tREF_NAME\t", r->name.c_str(), flag);
		// POS (1-based), MAPQ
		fprintf(samFile, "%d\t%d\t", (int)(r->aln.read_start+1), r->aln.score);

		// CIGAR
		fprintf(samFile, "CIGAR");

		// RNEXT, PNEXT, TLEN (print void mate position and coordinate)
		fprintf(samFile, "\t*\t0\t0\t");

		// SEQ, QUAL (print sequence and quality)
		const char* seq = r->aln.rc ? r->rc.c_str() : r->seq.c_str();
		for (int i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", "AGCTN"[(int)seq[i]]);
		}
		fprintf(samFile, "\t");
		for (int i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", '*');
		}
		fprintf(samFile, "\n");
	} else { // unmapped read
		int flag = SAM_FSU;
		// QNAME, FLAG
		fprintf(samFile, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", r->name.c_str(), flag);

		// SEQ, QUAL (print sequence and quality)
		for (int i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", "AGCTN"[(int)r->seq[i]]);
		}
		fprintf(samFile, "\t");
		for (int i = 0; i != r->len; i++) {
			fprintf(samFile, "%c", '*');
		}
		fprintf(samFile, "\n");
	}
}
