#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_matrix.h>

#include <netcdf.h>

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
ncdf_write_file(const char *filename, const struct gibbs_problem *p)
{
    int ncid;
    int time_dim_id;
    int cons_dim_id;
    int draw_dim_id;
    int dim_ids[3];
    int data_var_id;
    int mthet_var_id;
    int msig_var_id;
    size_t i;

    nc_create(filename, NC_CLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL, &ncid);

    nc_def_dim(ncid, "time", p->data->size1, &time_dim_id);
    nc_def_dim(ncid, "constituents", p->data->size2, &cons_dim_id);
    nc_def_dim(ncid, "draws", p->draws, &draw_dim_id);

    nc_put_att_long(ncid, NC_GLOBAL, "iterations", NC_INT, 1, (long *) &p->iterations);
    nc_put_att_long(ncid, NC_GLOBAL, "burn", NC_INT, 1, (long *) &p->burn);
    nc_put_att_long(ncid, NC_GLOBAL, "seed", NC_INT, 1, &p->seed);

    nc_def_var(ncid, "mean", NC_DOUBLE, 1, &cons_dim_id, &mthet_var_id);

    dim_ids[0] = cons_dim_id;
    dim_ids[1] = cons_dim_id;
    nc_def_var(ncid, "covariance", NC_DOUBLE, 2, dim_ids, &msig_var_id);

    dim_ids[2] = cons_dim_id;
    dim_ids[1] = time_dim_id;
    dim_ids[0] = draw_dim_id;
    nc_def_var(ncid, "data", NC_DOUBLE, 3, dim_ids, &data_var_id);
    nc_def_var_deflate(ncid, data_var_id, 1, 1, 5);

    nc_enddef(ncid);

    nc_put_var_double(ncid, mthet_var_id, p->mthet->data);

    if (p->msig->size2 == p->msig->tda) {
        nc_put_var_double(ncid, msig_var_id, p->msig->data);
    } else {
        for (i = 0; i < p->msig->size1; i++) {
            gsl_vector_const_view v = gsl_matrix_const_row(p->msig, i);
            size_t start[2] = {0, i};
            size_t count[2] = {p->msig->size2, 1};

            nc_put_vara_double(ncid, msig_var_id, start, count, v.vector.data);
        }
    }

    for (i = 0; i < p->draws; i++) {
        const gsl_matrix *m = p->ddata[i];

        if (m->size2 == m->tda) {
            size_t start[3] = {i, 0, 0};
            size_t count[3] = {1, m->size1, m->size2};

            nc_put_vara_double(ncid, data_var_id, start, count, m->data);
        } else {
            size_t j;

            for (j = 0; j < m->size1; i++) {
                gsl_vector_const_view v = gsl_matrix_const_row(m, j);
                size_t start[3] = {0, j, i};
                size_t count[3] = {m->size2, 1, 1};

                nc_put_vara_double(ncid, data_var_id, start, count, v.vector.data);
            }
        }
    }

    nc_close(ncid);
}

static void
print_usage(const char *progname)
{
    fprintf(stderr,
            "Usage: %s [OPTION] DATA MDLS OUT\n"
            "  -n int   number of iterations (default 1000)\n"
            "  -b int   number of iterations to burn-in (default 0)\n"
            "  -d int   number of random draws (default 1)\n"
            "  -r int   random seed (default 0)\n"
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
    const char *output_filename;
    struct gibbs_problem p = {
        NULL,
        NULL,
        1000,
        0,
        1,
        0,
        NULL, NULL, NULL, NULL, NULL};

    gsl_matrix *data = NULL;
    gsl_matrix *mdls = NULL;

    int optflag;
    int exit_val = 0;

    opterr = 0;

    while ((optflag = getopt(argc, argv, "n:b:d:r:h")) != -1) {
        char *endptr = NULL;

        switch (optflag) {
            case 'n':
                p.iterations = strtol(optarg, &endptr, 10);
                break;
            case 'b':
                p.burn = strtol(optarg, &endptr, 10);
                break;
            case 'd':
                p.draws = strtol(optarg, &endptr, 10);
                break;
            case 'r':
                p.seed = strtol(optarg, &endptr, 0);
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
        fprintf(stderr, "%s: missing argument OUT\n", argv[0]);
        print_usage(argv[0]);
        return 1;
    }
    output_filename = argv[optind++];

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

    p.data = data;
    p.mdls = mdls;

    gibbs_problem_exec(&p);
    ncdf_write_file(output_filename, &p);
    gibbs_problem_free(&p);

exit:
    if (data != NULL) {
        gsl_matrix_free(data);
    }

    if (mdls != NULL) {
        gsl_matrix_free(mdls);
    }

    return exit_val;
}
