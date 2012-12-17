#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stddef.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

void vector_remove_elem(const gsl_vector *u, size_t i, gsl_vector *v);
void matrix_row_remove_elem(const gsl_matrix *A, size_t i, size_t j, gsl_vector *v);
void matrix_remove_rowcol(const gsl_matrix *A, size_t r1, size_t r2, gsl_matrix *B);
void matrix_update_const_diag(gsl_matrix *A, double alpha, double beta);

__END_DECLS

#endif
