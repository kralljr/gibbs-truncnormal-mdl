#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include "csv.h"
#include "wishart.h"

static void
print_usage(const char *progname)
{
    fprintf(stderr,
            "Usage: %s [OPTION]\n"
            "  -n int               number of draws (default 1)\n"
            "  -d float             degrees of freedom (default 1.0)\n"
            "  -s float,float,...   scale\n"
            "  -o filename          output file (default standard output)\n"
            "  -r int               random seed (default 0)\n"
            "  -i                   draw from inverse Wishart\n"
            "  -h                   show this help\n"
            , progname);
}

static void
print_matrix(FILE *f, const gsl_matrix *m)
{
    size_t i;
    size_t j;

    fprintf(f, "%0.16f", gsl_matrix_get(m, 0, 0));
    for (j = 1; j < m->size2; j++) {
        fprintf(f, ",%0.16f", gsl_matrix_get(m, 0, j));
    }

    for (i = 1; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
            fprintf(f, ",%0.16f", gsl_matrix_get(m, i, j));
        }
    }
    fprintf(f, "\n");
}

int
main(int argc, char *argv[])
{
    int c;
    int n = 1;
    double d = 1.0;
    const char *scale_opt = NULL;
    double *scale_array;
    gsl_matrix_view scale;
    gsl_matrix *tmp_scale;
    gsl_matrix *z;
    size_t size;
    const char *filename = NULL;
    FILE *file;
    long seed = 0;
    int inverse = 0;
    int i;
    gsl_rng *rng;

    opterr = 0;

    while ((c = getopt(argc, argv, "n:d:s:o:r:ih")) != -1) {
        char *endptr = NULL;

        switch (c) {
            case 'n':
                n = strtol(optarg, &endptr, 10);
                break;
            case 'd':
                d = strtod(optarg, &endptr);
                break;
            case 's':
                scale_opt = optarg;
                break;
            case 'o':
                filename = optarg;
                break;
            case 'r':
                seed = strtol(optarg, &endptr, 0);
                break;
            case 'i':
                inverse = 1;
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

    if (scale_opt == NULL) {
        fprintf(stderr, "%s: missing scale\n", argv[0]);
        print_usage(argv[0]);
        return 3;
    }

    size = csv_count_fields(scale_opt);
    size = (int) sqrt((double) size);

    scale_array = calloc(size * size, sizeof (double));
    csv_read_double_fields(scale_array, size * size, scale_opt);
    scale = gsl_matrix_view_array(scale_array, size, size);

    tmp_scale = gsl_matrix_alloc(size, size);
    z = gsl_matrix_alloc(size, size);

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

    for (i = 0; i < n; i++) {
        gsl_matrix_memcpy(tmp_scale, &scale.matrix);
        if (inverse) {
            ran_invwishart(rng, d, tmp_scale, z);
        } else {
            ran_wishart(rng, d, tmp_scale, z);
        }
        print_matrix(file, z);
    }

    if (filename != NULL) {
        fclose(file);
    }

    free(scale_array);
    gsl_matrix_free(tmp_scale);
    gsl_matrix_free(z);
    gsl_rng_free(rng);

    return 0;
}
