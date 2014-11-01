CC=gcc
NC_CFLAGS=$(shell pkg-config --cflags netcdf)
NC_LDFLAGS=$(shell pkg-config --libs netcdf)
GSL_CFLAGS=$(shell pkg-config --cflags gsl)
GSL_LDFLAGS=$(shell pkg-config --libs gsl)
CFLAGS=$(NC_CFLAGS) $(GSL_CFLAGS) -Wall -Wextra -O3 -DHAVE_INLINE
LDFLAGS=$(NC_LDFLAGS) $(GSL_LDFLAGS) -lm

# If using atlas, uncomment the following line to enable
#LDFLAGS=$(NC_LDFLAGS) -lgsl -L/usr/lib64/atlas -lcblas -latlas -lm

# On OS X, uncomment the following line to enable the Accelerate framework
#LDFLAGS=$(NC_LDFLAGS) -lgsl -framework Accelerate -lm

# Version
VERSION = 1.0


all: gibbs rwishart rmvnorm rtnorm

gibbs: csv.o gibbs.o main.o matrix.o multivariate.o mvnormal.o truncnormal.o wishart.o
	$(CC) $(LDFLAGS) $^ -o $@

rwishart: rwishart.c csv.o wishart.o
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


# Phony targets

dist:
	mkdir gibbs-truncnormal-mdl-$(VERSION)
	cp *.c *.h Makefile README.md gibbs-truncnormal-mdl-$(VERSION)
	tar czf gibbs-truncnormal-mdl-$(VERSION).tar.gz gibbs-truncnormal-mdl-$(VERSION)

clean:
	-rm -rf gibbs-truncnormal-mdl-$(VERSION)
	-rm -f *.o gibbs rwishart rmvnorm rtnorm
