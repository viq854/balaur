CC=		g++
CFLAGS=		-g -Wall -O3 -fopenmp -std=c++11
DFLAGS=		-DHAVE_PTHREAD
OBJS=		main.o io.o index.o hash.o city.o cluster.o align.o mt19937-64.o
PROG=		srx
INCLUDES=	
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o .cc

.cc.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)
		
srx:$(OBJS)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
