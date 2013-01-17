#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_matrix.h>

#include "csv.h"
#include "gibbs.h"

/* XXX: use glibc's getline if available */
static ssize_t
fgets_checked(char *s, size_t n, FILE *f)
{
    size_t len;

    if (n > 0) {
        *s = '\0';
    }

    fgets(s, n, f);
    len = strlen(s);

    if (len > 1 && s[len - 1] != '\n') {
        return -1;
    }

    return len;
}
 
/* XXX: This whole 2-pass approach is a disaster */
static gsl_matrix *
csv_read_file(FILE *f)
{
    char *lineptr;
    size_t linelen = 16384;
    size_t nfields;
    size_t nlines;
    long pos;
    double *fields;
    gsl_matrix *matrix;
    size_t i;

    lineptr = malloc(linelen);

    pos = ftell(f);
    nlines = 0;
    while (!ferror(f) && !feof(f)) {
        ssize_t rval = fgets_checked(lineptr, linelen, f);

        if (rval == -1) {
            fprintf(stderr, "Line too long on line: %d\n", (int) nlines);
        }

        for (i = 0; i < strlen(lineptr); i++) {
            if (!isspace(lineptr[i])) {
                nlines++;
                break;
            }
        }
    }
    fseek(f, pos, SEEK_SET);

    /* Assume that the first line is a header */
    fgets(lineptr, linelen, f);
    nfields = csv_count_fields(lineptr);
    fields = calloc(nfields, sizeof(double));

    matrix = gsl_matrix_alloc(nlines - 1, nfields);
    for (i = 0; i < nlines - 1; i++) {
        gsl_vector_view row;
        gsl_vector_view in;
        const char *parseval;

        fgets(lineptr, linelen, f);
        parseval = csv_read_double_fields(fields, nfields, lineptr);
        if (parseval != NULL) {
            fprintf(stderr, "Unexpected character on line %d: at position: %d",
                    (int) i + 1, (int) (parseval - lineptr));
            return NULL;
        }

        row = gsl_matrix_row(matrix, i);
        in = gsl_vector_view_array(fields, nfields);
        gsl_vector_memcpy(&row.vector, &in.vector);
    }

    free(lineptr);
    free(fields);

    return matrix;
}

static gsl_matrix *
csv_load_file(const char *filename)
{
    FILE *f;
    gsl_matrix *M;

    f = fopen(filename, "rt");
    if (f == NULL) {
        return NULL;
    }

    M = csv_read_file(f);
    fclose(f);

    return M;
}

static void
print_usage(const char *progname)
{
    fprintf(stderr,
            "Usage: %s [OPTION] DATA MDLS OUTDIR\n"
            "  -n int   number of iterations (default 1000)\n"
            "  -b int   number of iterations to burn-in (default 0)\n"
            "  -d int   number of random draws (default 1)\n"
            "  -r int   random seed (default 0)\n"
            "  -p int   print progress every n iterations (default 1000)\n"
            "  -h       show this help\n"
            , progname);
}

/* TODO: argument parsing is rudimentary.  Need to validate numerical
 * arguments. */
int
main(int argc, char * const argv[])
{
    const char *data_filename;
    const char *mdls_filename;
    const char *output_directory;
    size_t iterations = 1000;
    size_t burn = 0;
    size_t draws = 1;
    size_t progress = 1000;
    long seed = 0;

    gsl_matrix *data = NULL;
    gsl_matrix *mdls = NULL;

    int optflag;
    int exit_val = 0;

    opterr = 0;

    while ((optflag = getopt(argc, argv, "n:b:d:r:p:h")) != -1) {
        char *endptr = NULL;

        switch (optflag) {
            case 'n':
                iterations = strtol(optarg, &endptr, 10);
                break;
            case 'b':
                burn = strtol(optarg, &endptr, 10);
                break;
            case 'd':
                draws = strtol(optarg, &endptr, 10);
                break;
            case 'r':
                seed = strtol(optarg, &endptr, 0);
                break;
            case 'p':
                progress = strtol(optarg, &endptr, 0);
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case '?':
                fprintf(stderr, "%s: invalid option or missing argument: -'%c'\n",
                        argv[0], optopt);
                print_usage(argv[0]);
                return 1;
        }
    }

    if (optind == argc) {
        fprintf(stderr, "%s: missing argument DATA\n", argv[0]);
        print_usage(argv[0]);
        return 1;
    }
    data_filename = argv[optind++];

    if (optind == argc) {
        fprintf(stderr, "%s: missing argument MDLS\n", argv[0]);
        print_usage(argv[0]);
        return 1;
    }
    mdls_filename = argv[optind++];

    if (optind == argc) {
        fprintf(stderr, "%s: missing argument OUTDIR\n", argv[0]);
        print_usage(argv[0]);
        return 1;
    }
    output_directory = argv[optind++];

    if (optind < argc) {
        fprintf(stderr, "%s: unexpected argument: %s\n", argv[0], argv[optind]);
        print_usage(argv[0]);
        return 1;
    }

    data = csv_load_file(data_filename);
    if (data == NULL) {
        fprintf(stderr, "%s: cannot read '%s': %s\n",
                argv[0], data_filename, strerror(errno));
        exit_val = errno;
        goto exit;
    }

    mdls = csv_load_file(mdls_filename);
    if (mdls == NULL) {
        fprintf(stderr, "%s: cannot read '%s': %s\n",
                argv[0], mdls_filename, strerror(errno));
        exit_val = errno;
        goto exit;
    }

    if (data->size1 != mdls->size1 || data->size2 != mdls->size2) {
        fprintf(stderr, "Non-matching `data' and `mdls' dimensions\n");
        exit_val = GSL_EBADLEN;
        goto exit;
    }

    impute_data(data, mdls, output_directory, iterations, burn, draws, progress,
            seed);

exit:
    if (data != NULL) {
        gsl_matrix_free(data);
    }

    if (mdls != NULL) {
        gsl_matrix_free(mdls);
    }

    return exit_val;
}
