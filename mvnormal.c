#include <stddef.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

#include "mvnormal.h"

void
ran_multivariate_normal_chol(const gsl_rng *rng, const gsl_vector *mean,
        const gsl_matrix *sigmachol, gsl_vector *z)
{
    size_t i;

    for (i = 0; i < z->size; i++) {
        double r;

        r = gsl_ran_gaussian(rng, 1.0);
        gsl_vector_set(z, i, r);
    }

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, sigmachol, z);
    gsl_vector_add(z, mean);
}

void
ran_multivariate_normal(const gsl_rng *rng, const gsl_vector *mean,
        gsl_matrix *sigma, gsl_vector *z)
{
    gsl_linalg_cholesky_decomp(sigma);
    ran_multivariate_normal_chol(rng, mean, sigma, z);
}
