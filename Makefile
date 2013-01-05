CC=gcc
CFLAGS=-Wall -Wextra -O0 -g -pg
LDFLAGS=-lm -lgslcblas -lgsl -g -pg

all: gibbs rtnorm

gibbs: gibbs.o main.o matrix.o multivariate.o truncnormal.o wishart.o
	$(CC) $(LDFLAGS) $^ -o $@

rtnorm: rtnorm.c truncnormal.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

gibbs.o: gibbs.c gibbs.h
	$(CC) $(CFLAGS) -c $<

main.o: main.c
	$(CC) $(CFLAGS) -c $<

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c $<

multivariate.o: multivariate.c multivariate.h
	$(CC) $(CFLAGS) -c $<

truncnormal.o: truncnormal.c truncnormal.h
	$(CC) $(CFLAGS) -c $<

wishart.o: wishart.c wishart.h
	$(CC) $(CFLAGS) -c $<

clean:
	-rm -f *.o gibbs
