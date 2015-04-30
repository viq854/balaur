# TBB library path
ENABLE_TBB=1
TBB_LIBRARY_DIR=
TBB_INCLUDE_DIR=

# MARISA trie library
ENABLE_MARISA=0
MARISA_LIBRARY_DIR=	#/usr/local/lib
MARISA_INCLUDE_DIR= 	#/usr/local/include/marisa

PROG=		totoro
CC=		g++
CFLAGS=		-Wall -Ofast -fopenmp -std=c++0x
DFLAGS=		-DUSE_TBB=$(ENABLE_TBB) -DUSE_MARISA=$(ENABLE_MARISA)

SOURCES=	main.cc \
		lsh.cc \
		io.cc \
		index.cc \
		city.cc \
		align.cc \
		mt19937-64.cc \
		sam.cc \

DEPS=		index.h align.h io.h city.h lsh.h sam.h			
OBJDIR=		obj
_OBJS=		$(SOURCES:.cc=.o)
OBJS=		$(patsubst %,$(OBJDIR)/%,$(_OBJS))
OBJDIR_GUARD = @mkdir -p $(OBJDIR)

INCLUDES=	-I$(TBB_INCLUDE_DIR) -I$(MARISA_INCLUDE_DIR)
LIBS=		-lm -lz -L$(TBB_LIBRARY_DIR) -ltbb -ltbbmalloc -L$(MARISA_LIBRARY_DIR) -lmarisa


$(OBJDIR)/%.o: %.cc $(DEPS)
		$(OBJDIR_GUARD)
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
                
all: $(PROG)
		
$(PROG): $(OBJS)
		@mkdir -p $(OBJDIR)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)

debug: CFLAGS += -g
debug: $(PROG)

clean:
		rm -f gmon.out $(OBJDIR)/*.o a.out $(PROG) *~ *.a
