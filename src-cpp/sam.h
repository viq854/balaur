#ifndef __SAM_H
#define __SAM_H

#include "index.h"

void print_aln2sam(FILE* samFile, read_t* r, const ref_t& ref);

struct sam_writer_t {
	FILE* samFile;

	void open_file(const std::string& fname) {
		std::string samFname(fname);
		samFname += std::string(".sam");
		FILE* samFile = (FILE*) fopen(samFname.c_str(), "w");
		if (samFile == NULL) {
			printf("sam_writer_t open: Cannot open SAM file: %s!\n", samFname.c_str());
			exit(1);
		}
	}

	void write_sam_batch(reads_t& reads, const ref_t& ref) {
		for (uint32 i = 0; i < reads.reads.size(); i++) {
			read_t* r = &reads.reads[i];
			print_aln2sam(samFile, r, ref);
		}
	}
	
	void close_file() {
		if(samFile != NULL)  fclose(samFile);
	}
};

void store_alns_sam(reads_t& reads, const ref_t& ref, const index_params_t* params);

#endif
