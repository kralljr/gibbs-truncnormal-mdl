#ifndef __MULTIVARIATE_H__
#define __MULTIVARIATE_H__

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

void multivariate_mean(const gsl_matrix *A, gsl_vector *mu);
void multivariate_covariance(const gsl_matrix *A, const gsl_vector *mu,
        gsl_matrix *S, gsl_vector *tmp);

__END_DECLS

#endif
