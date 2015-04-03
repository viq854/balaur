#ifndef __SAM_H
#define __SAM_H

#include "index.h"

int openOutFile(char * filename);
int writeAlignment(read_t& r);
void closeOutFile();

#endif
