#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "wishart.h"

void
ran_wishart(const gsl_rng *r, double nu, gsl_matrix *scale, gsl_matrix *rwish)
{
    size_t i;
    size_t j;

    /* TODO: check that nu >= scale->size */

    /* A <- AA' */
    gsl_linalg_cholesky_decomp(scale);

    /* B_{ij} = 0 */
    gsl_matrix_set_zero(rwish);

    /* B_{ii} ~ Gamma(k - i + 1, 1/2), 1 <= i <= p */
    for (i = 1; i <= rwish->size1; i++) {
        double y;

        y = sqrt(gsl_ran_gamma(r, (nu - i + 1) / 2.0, 2.0));
        gsl_matrix_set(rwish, i - 1, i - 1, y);
    }

    /* B_{ij} ~ N(0, 1), 1 <= j < i <= p */
    for (i = 1; i <= rwish->size1; i++) {
        for (j = 1; j < i; j++) {
            double y;

            y = gsl_ran_ugaussian(r);
            gsl_matrix_set(rwish, i - 1, j - 1, y);
        }
    }

    /* B <- AB */
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0,
            scale, rwish);

    /* A <- BB' */
    gsl_matrix_memcpy(scale, rwish);
    gsl_blas_dtrmm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1.0,
            scale, rwish);
}

void
ran_invwishart(const gsl_rng *r, double nu, gsl_matrix *scale,
        gsl_matrix *riwish)
{
    gsl_linalg_cholesky_decomp(scale);
    gsl_linalg_cholesky_invert(scale);

    ran_wishart(r, nu, scale, riwish);

    gsl_linalg_cholesky_decomp(riwish);
    gsl_linalg_cholesky_invert(riwish);
}
