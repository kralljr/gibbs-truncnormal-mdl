#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "wishart.h"

void
ran_wishart(const gsl_rng *r, double nu, gsl_matrix *chol, gsl_matrix *tmp)
{
    size_t i;
    size_t j;

    /* A <- AA' */
    gsl_linalg_cholesky_decomp(chol);

    /* B_{ij} = 0 */
    gsl_matrix_set_zero(tmp);

    /* B_{ii} ~ Gamma(k - i + 1, 1/2), 1 <= i <= p */
    for (i = 1; i <= tmp->size1; i++) {
        double y;

        y = gsl_ran_gamma(r, nu - i + 1, 0.5);
        gsl_matrix_set(tmp, i - 1, i - 1, y);
    }

    /* B_{ij} ~ N(0, 1), 1 <= j < i <= p */
    for (i = 1; i <= tmp->size1; i++) {
        for (j = 1; j < i; j++) {
            double y;

            y = gsl_ran_ugaussian(r);
            gsl_matrix_set(tmp, i - 1, j - 1, y);
        }
    }

    /* B <- AB */
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, chol,
            tmp);

    /* A <- BB' */
    gsl_matrix_memcpy(chol, tmp);
    gsl_blas_dtrmm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1.0, chol,
            tmp);
}

void
ran_invwishart(const gsl_rng *r, double nu, gsl_matrix *chol,
        gsl_permutation *permutation, gsl_matrix *riwish)
{
    int signum;

    gsl_linalg_cholesky_decomp(chol);
    gsl_linalg_cholesky_invert(chol);
    ran_wishart(r, nu, chol, riwish);

    gsl_linalg_LU_decomp(riwish, permutation, &signum);
    gsl_matrix_memcpy(chol, riwish);
    gsl_linalg_LU_invert(chol, permutation, riwish);
}
