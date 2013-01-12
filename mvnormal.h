#ifndef __MVNORMAL_H__
#define __MVNORMAL_H__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
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

void ran_multivariate_normal_chol(const gsl_rng *rng, const gsl_vector *mean,
        const gsl_matrix *sigmachol, gsl_vector *z);
void ran_multivariate_normal(const gsl_rng *rng, const gsl_vector *mean,
        gsl_matrix *sigma, gsl_vector *z);

__END_DECLS

#endif
