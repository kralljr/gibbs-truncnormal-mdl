CC=gcc
CFLAGS=-Wall -Wextra -O2 -DHAVE_INLINE
LDFLAGS=-lm -lgslcblas -lgsl

all: gibbs riwishart rmvnorm rtnorm

gibbs: csv.o gibbs.o main.o matrix.o multivariate.o mvnormal.o truncnormal.o wishart.o
	$(CC) $(LDFLAGS) $^ -o $@

riwishart: riwishart.c csv.o wishart.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

rmvnorm: rmvnorm.c csv.o mvnormal.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

rtnorm: rtnorm.c truncnormal.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

csv.o: csv.c csv.h
	$(CC) $(CFLAGS) -c $<

gibbs.o: gibbs.c gibbs.h
	$(CC) $(CFLAGS) -c $<

main.o: main.c
	$(CC) $(CFLAGS) -c $<

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c $<

multivariate.o: multivariate.c multivariate.h
	$(CC) $(CFLAGS) -c $<

mvnormal.o: mvnormal.c mvnormal.h
	$(CC) $(CFLAGS) -c $<

truncnormal.o: truncnormal.c truncnormal.h
	$(CC) $(CFLAGS) -c $<

wishart.o: wishart.c wishart.h
	$(CC) $(CFLAGS) -c $<

clean:
	-rm -f *.o gibbs
