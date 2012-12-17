CC=gcc
CFLAGS=-Wall -Wextra -O0 -g
LDFLAGS=-lm -lgslcblas -lgsl

all: gibbs

gibbs: gibbs.o main.o matrix.o multivariate.o wishart.o
	$(CC) $(LDFLAGS) $^ -o $@

gibbs.o: gibbs.c gibbs.h
	$(CC) $(CFLAGS) -c $<

main.o: main.c
	$(CC) $(CFLAGS) -c $<

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c $<

multivariate.o: multivariate.c multivariate.h
	$(CC) $(CFLAGS) -c $<

wishart.o: wishart.c wishart.h
	$(CC) $(CFLAGS) -c $<

clean:
	-rm -f *.o gibbs
