#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_matrix.h>

#include "csv.h"
#include "gibbs.h"

/* FIXME: use glibc's getline if available */
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
        perror(filename);
 
        return NULL;
    }

    M = csv_read_file(f);
    fclose(f);

    return M;
}

int
main(int argc, const char *argv[])
{
    gsl_matrix *data;
    gsl_matrix *mdls;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s DATA MDLS\n", argv[0]);
        
        return EINVAL;
    }

    data = csv_load_file(argv[1]);
    if (data == NULL) {
        return errno;
    }

    mdls = csv_load_file(argv[2]);
    if (data == NULL) {
        return errno;
    }

    if (data->size1 != mdls->size1 || data->size2 != mdls->size2) {
        fprintf(stderr, "Non-matching `data' and `mdls' dimensions\n");

        return GSL_EBADLEN;
    }

    impute_data(data, mdls);

    gsl_matrix_free(data);
    gsl_matrix_free(mdls);

    return 0;
}
