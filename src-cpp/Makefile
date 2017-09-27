# Dependencies

# OpenMP currently required (Seqan compatibility issue) 
ENABLE_OPENMP=1

# TBB library path (recommended)
ENABLE_TBB=1
TBB_INCLUDE_DIR=
TBB_LIBRARY_DIR=

# Seqan (Fastq IO etc)
SEQAN_INCLUDE_DIR=../third-party/seqan/include/
SEQAN_OPENMP=$(ENABLE_OPENMP)
SEQAN_ZLIB=1

# Encryption
USE_SHA1_ASM=1
ASM_SOURCES=../third-party/sha1-fast-64.S

# Development-only experimental/temporary libraries
ENABLE_MARISA=0
#MARISA_LIBRARY_DIR=#/usr/local/lib
#MARISA_INCLUDE_DIR=#/usr/local/include/marisa 
#COMP_LIB=/home/viq/libvbyte
#OROCH_LIB=/home/viq/Oroch/oroch
#LEMON_INCLUDE_DIR=/home/viq/lemon/include
#LEMON_LIB_DIR=/home/viq/lemon/lib
#CLASP_LIB=/home/viq/balaur/src-cpp/clasp_v1_1/libs
#LOCAL_INCLUDE_DIR=/home/viq/local/include
#LOCAL_LIB_DIR=/home/viq/local/lib

PROG=		balaur
CC=		g++
CFLAGS=		-g -msse4.1 -o rand_sse -O3 -std=c++0x
DFLAGS=		-DSEQAN_HAS_ZLIB=$(SEQAN_ZLIB) -DSEQAN_HAS_OPENMP=$(SEQAN_OPENMP) -DUSE_SHA1_ASM=$(USE_SHA1_ASM) -DUSE_TBB=$(ENABLE_TBB) -DUSE_MARISA=$(ENABLE_MARISA)

ifeq ($(ENABLE_OPENMP), 1)
CFLAGS   +=-fopenmp
endif

SOURCES=main.cc \
	lsh.cc \
	io.cc \
	index.cc \
	align.cc \
	sam.cc \
	stats.cc \
	crypt.cc \
	voting.cc \
	contigs.cc \
	eval.cc \
	mt19937-64.cc \
	city.cc	

DEPS=		index.h align.h io.h lsh.h sam.h hash.h voting.h contigs.h crypt.h			
OBJDIR=		obj
_OBJS=		$(SOURCES:.cc=.o)
OBJS=		$(patsubst %,$(OBJDIR)/%,$(_OBJS))
OBJDIR_GUARD = @mkdir -p $(OBJDIR)

INCLUDES=	-I/usr/local/include -I$(SEQAN_INCLUDE_DIR)
ifneq ($(TBB_INCLUDE_DIR),)
INCLUDES+= -I$(TBB_INCLUDE_DIR)
endif
LIBS=		-L/usr/local/lib -lz -lm -lrt
ifneq ($(TBB_LIBRARY_DIR),)
LIBS+= -L$(TBB_LIBRARY_DIR)
endif
ifeq ($(ENABLE_TBB), 1)
LIBS+= -ltbb -ltbbmalloc
endif

#INCLUDES_EXP=       -I$(COMP_LIB) -I$(LOCAL_INCLUDE_DIR)-I$(MARISA_INCLUDE_DIR) -I$(LEMON_INCLUDE_DIR)
#LIBS_EXP=           -lcrypto -lm -lz -L$(MARISA_LIBRARY_DIR) -lmarisa -lrt -L$(CLASP_LIB) -lclasp -L${LEMON_LIB_DIR} -lemon -L$(LOCAL_LIB_DIR) -lb2 -L$(COMP_LIB) -lvbyte

$(OBJDIR)/%.o: %.cc $(DEPS)
		$(OBJDIR_GUARD)
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
                
all: $(PROG)
		
$(PROG): $(OBJS)
		@mkdir -p $(OBJDIR)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) $(ASM_SOURCES) -o $@ $(LIBS)

debug: CFLAGS += -g -DDEBUG
debug: $(PROG)

sim: CFLAGS += -DSIM_EVAL
sim: $(PROG)

vanilla: CFLAGS += -DSIM_EVAL -DVANILLA
vanilla: $(PROG)

clean:
		rm -f gmon.out $(OBJDIR)/*.o a.out $(PROG) *~ *.a
