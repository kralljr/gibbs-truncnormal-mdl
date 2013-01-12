#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

#include "csv.h"
#include "mvnormal.h"

static void
print_usage(const char *progname)
{
    fprintf(stderr,
            "Usage: %s [OPTION]\n"
            "  -n int               number of draws (default 1)\n"
            "  -m float,float,...   mean\n"
            "  -c float,float,...   covariance\n"
            "  -o filename          output file (default standard output)\n"
            "  -r int               random seed (default 0)\n"
            "  -h                   show this help\n"
            , progname);
}

static void
print_vector(FILE *f, const gsl_vector *v)
{
    size_t i;

    fprintf(f, "%0.16f", gsl_vector_get(v, 0));
    for (i = 1; i < v->size; i++) {
        fprintf(f, ",%0.16f", gsl_vector_get(v, i));
    }
    fprintf(f, "\n");
}

int
main(int argc, char *argv[])
{
    int c;
    int n = 1;
    const char *mean_opt = NULL;
    const char *covariance_opt = NULL;
    double *mean_array;
    double *covariance_array;
    gsl_vector_view mean;
    gsl_matrix_view covariance;
    gsl_vector *z;
    size_t size;
    const char *filename = NULL;
    FILE *file;
    long seed = 0;
    int i;
    gsl_rng *rng;

    opterr = 0;

    while ((c = getopt(argc, argv, "n:m:c:o:r:h")) != -1) {
        char *endptr = NULL;

        switch (c) {
            case 'n':
                n = strtol(optarg, &endptr, 10);
                break;
            case 'm':
                mean_opt = optarg;
                break;
            case 'c':
                covariance_opt = optarg;
                break;
            case 'o':
                filename = optarg;
                break;
            case 'r':
                seed = strtol(optarg, &endptr, 0);
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case '?':
                fprintf(stderr, "%s: invalid option or missing argument: -%c\n",
                        argv[0], optopt);
                print_usage(argv[0]);
                return 1;
        }
    }

    if (mean_opt == NULL) {
        fprintf(stderr, "%s: missing mean\n", argv[0]);
        print_usage(argv[0]);
        return 3;
    }

    if (covariance_opt == NULL) {
        fprintf(stderr, "%s: missing covariance\n", argv[0]);
        print_usage(argv[0]);
        return 3;
    }

    size = csv_count_fields(mean_opt);

    mean_array = calloc(size, sizeof (double));
    csv_read_double_fields(mean_array, size, mean_opt);
    mean = gsl_vector_view_array(mean_array, size);

    covariance_array = calloc(size * size, sizeof (double));
    csv_read_double_fields(covariance_array, size * size, covariance_opt);
    covariance = gsl_matrix_view_array(covariance_array, size, size);

    z = gsl_vector_alloc(size);

    if (filename == NULL) {
        file = stdout;
    } else {
        file = fopen(filename, "wt");
        if (file == NULL) {
            fprintf(stderr, "%s: can't open %s: %s\n", argv[0], filename,
                    strerror(errno));
            return 2;
        }
    }

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, seed);

    ran_multivariate_normal(rng, &mean.vector, &covariance.matrix, z);
    print_vector(file, z);
    for (i = 1; i < n; i++) {
        ran_multivariate_normal_chol(rng, &mean.vector, &covariance.matrix, z);
        print_vector(file, z);
    }

    if (filename != NULL) {
        fclose(file);
    }

    free(mean_array);
    free(covariance_array);
    gsl_vector_free(z);
    gsl_rng_free(rng);

    return 0;
}
