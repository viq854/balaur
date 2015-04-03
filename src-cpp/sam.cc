#include "sam.h"

FILE * outfile = NULL;

int openOutFile(char * filename)
{
  outfile = fopen(filename, "w");

  return 0;
}

int writeAlignment(read_t& r)
{
  if(outfile)
  {
    fprintf(outfile, "%s\t0\tREF_NAME\t%d\t255\t---CIGAR STRING---\tRNEXT\t\
PNEXT\tTLEN\t%s\t%s\n", r.name.c_str(), r.ref_pos_l, r.seq.c_str(),
    r.seq.c_str());
  }

  return 0;
}

void closeOutFile()
{
  if (outfile)
    fclose(outfile);
}
