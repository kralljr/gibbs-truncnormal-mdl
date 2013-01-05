#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_rng.h>

#include "truncnormal.h"

static void
exit_usage(char *argv[])
{
    fprintf(stderr,
            "Usage: %s [OPTION]\n"
            "  -n int       number of draws (default 1)\n"
            "  -m float     mean (default 0.0)\n"
            "  -s float     standard deviation (default 1.0)\n"
            "  -a float     lower bound (default -1.0)\n"
            "  -b float     upper bound (default 1.0)\n"
            "  -o filename  output file (default standard output)\n"
            "  -r int       random seed (default 0)\n"
            "  -h           show this help\n"
            , argv[0]);

    exit(0);
}

static void
exit_eio(char *argv[], char *filename)
{
    fprintf(stderr, "%s: can't open %s: %s\n", argv[0], filename, strerror(errno));
    exit(2);
}

static void
exit_invalid_flag(char *argv[], char flag)
{
    fprintf(stderr, "%s: invalid option or missing argument: -%c\n", argv[0], flag);
    exit(1);
}

static void
exit_invalid_argument(char *argv[], char flag, char *arg)
{
    fprintf(stderr, "%s: invalid argument to option -%c: %s\n", argv[0], flag, arg);
    exit(1);
}

static double
checked_strtod(char *argv[], char flag, char *s)
{
    double y;
    char *endptr = NULL;

    y = strtod(s, &endptr);
    if (s == endptr) {
        exit_invalid_argument(argv, flag, s);
    }

    return y;
}

int main(int argc, char *argv[])
{
    int c;
    int n = 1;
    double m = 0.0;
    double s = 1.0;
    double a = -1.0;
    double b = 1.0;
    char *filename = NULL;
    FILE *file;
    long seed = 0;
    int i;
    gsl_rng *rng;

    opterr = 0;

    while ((c = getopt(argc, argv, "n:m:s:a:b:o:r:h")) != -1) {
        char *endptr = NULL;

        switch (c) {
            case 'n':
                n = strtol(optarg, &endptr, 10);
                break;
            case 'm':
                m = checked_strtod(argv, 'm', optarg);
                break;
            case 's':
                s = checked_strtod(argv, 's', optarg);
                break;
            case 'a':
                a = checked_strtod(argv, 'a', optarg);
                break;
            case 'b':
                b = checked_strtod(argv, 'b', optarg);
                break;
            case 'o':
                filename = optarg;
                break;
            case 'r':
                seed = strtol(optarg, &endptr, 0);
                break;
            case 'h':
                exit_usage(argv);
                break;
            case '?':
                exit_invalid_flag(argv, optopt);
                break;
        }
    }

    if (filename == NULL) {
        file = stdout;
    } else {
        file = fopen(filename, "wt");
        if (file == NULL) {
            exit_eio(argv, filename);
        }
    }

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, seed);

    for (i = 0; i < n; i++) {
        double y;

        y = ran_truncnormal(rng, a, b, m, s);

        fprintf(file, "%0.16f\n", y);
    }

    if (filename != NULL) {
        fclose(file);
    }

    gsl_rng_free(rng);

    return 0;
}
