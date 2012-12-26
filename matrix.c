/**
 * Additional matrix operations.
 */

#include <stddef.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "matrix.h"

/**
 * Subvector with the r-th element removed, i.e., v <- u[-i].
 */
void
vector_remove_elem(const gsl_vector *u, size_t i, gsl_vector *v)
{
    if (i > 0) {
        gsl_vector_const_view uhead = gsl_vector_const_subvector(u, 0, i);
        gsl_vector_view vhead = gsl_vector_subvector(v, 0, i);
        gsl_vector_memcpy(&vhead.vector, &uhead.vector);
    }

    if (i + 1 < u->size) {
        gsl_vector_const_view utail =
                gsl_vector_const_subvector(u, i + 1, u->size - i - 1);
        gsl_vector_view vtail =
                gsl_vector_subvector(v, i, u->size - i - 1);
        gsl_vector_memcpy(&vtail.vector, &utail.vector);
    }
}

/**
 * The i'th row of a matrix with the j'th element removed, i.e., v <- A[i, -j].
 */
void
matrix_row_remove_elem(const gsl_matrix *A, size_t i, size_t j, gsl_vector *v)
{
    gsl_vector_const_view u = gsl_matrix_const_row(A, i);
    vector_remove_elem(&u.vector, j, v);
}

/**
 * Submatrix with the r1-th row and r2-th column removed, i.e., B <- A[-i, -j].
 */
void
matrix_remove_rowcol(const gsl_matrix *A, size_t r1, size_t r2, gsl_matrix *B)
{
    size_t i;
    size_t j = 0;

    for (i = 0; i < r1; i++, j++) {
        gsl_vector_view v = gsl_matrix_row(B, j);
        matrix_row_remove_elem(A, i, r2, &v.vector);
    }

    for (i = r1 + 1; i < A->size1; i++, j++) {
        gsl_vector_view v = gsl_matrix_row(B, j);
        matrix_row_remove_elem(A, i, r2, &v.vector);
    }
}

/**
 * For scalars alpha and beta, A <- alpha * A + beta * I
 */
void
matrix_update_const_diag(gsl_matrix *A, double alpha, double beta)
{
    size_t i;

    if (alpha == 1.0) {
        for (i = 0; i < A->size1; i++) {
            double y = gsl_matrix_get(A, i, i);

            gsl_matrix_set(A, i, i, y + beta);
        }
    } else if (alpha == -1.0) {
        for (i = 0; i < A->size1; i++) {
            double y = gsl_matrix_get(A, i, i);

            gsl_matrix_set(A, i, i, beta - y);
        }
    } else {
        for (i = 0; i < A->size1; i++) {
            double y = gsl_matrix_get(A, i, i);

            gsl_matrix_set(A, i, i, alpha * y + beta);
        }
    }
}

/**
 * Compute the inverse of A with a row and column removed, given the inverse
 * of A.  If A is n-by-n, then B is (n-1)-by-(n-1), and work is (n-1).
 *
 * Computes inv(A[-i, -i]) given inv(A) for symmetric A.
 */
void
matrix_invert_remove_rowcol(const gsl_matrix *Ainv, size_t i, gsl_matrix *Binv,
        gsl_vector *work)
{
    matrix_remove_rowcol(Ainv, i, i, Binv);
    matrix_row_remove_elem(Ainv, i, i, work);
    gsl_blas_dger(-1.0 / gsl_matrix_get(Ainv, i, i), work, work, Binv);
}
