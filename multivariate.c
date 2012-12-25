#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "multivariate.h"

void
multivariate_mean(const gsl_matrix *A, gsl_vector *mu)
{
    size_t i;

    gsl_vector_set_zero(mu);
    for (i = 0; i < A->size1; i++) {
        gsl_vector_const_view row = gsl_matrix_const_row(A, i);
        gsl_vector_add(mu, &row.vector);
    }

    gsl_vector_scale(mu, 1.0 / A->size1);
}

void
multivariate_covariance(const gsl_matrix *A, const gsl_vector *mu,
        gsl_matrix *S, gsl_vector *tmp)
{
    size_t i;

    gsl_matrix_set_zero(S);
    for (i = 0; i < A->size1; i++) {
        gsl_vector_const_view row = gsl_matrix_const_row(A, i);
        gsl_vector_memcpy(tmp, &row.vector);
        gsl_vector_sub(tmp, mu);

        gsl_blas_dger(1.0 / A->size1, tmp, tmp, S);
    }
}
