CC=		g++
CFLAGS=		-g -Wall -O3 -fopenmp -std=c++11
DFLAGS=		-DHAVE_PTHREAD
OBJS=		main.o io.o index.o lsh.o city.o cluster.o align.o mt19937-64.o ksw.o
PROG=		srx
INCLUDES=	-I/usr/local/include/ -I/home/viq/seqan/core/include/ -I/home/viq/seqan/extras/include/
LIBS=		-lm -lz -ltbb -ltbbmalloc -lpthread -L/usr/local/bin -lmarisa -lrt

.SUFFIXES:.c .o .cc

.cc.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

.c.o:
		gcc -c -g  $(INCLUDES) $< -o $@

all:$(PROG)
		
srx:$(OBJS)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
