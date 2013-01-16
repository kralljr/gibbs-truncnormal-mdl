#ifndef __WISHART_H__
#define __WISHART_H__

#include <stddef.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

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

void ran_wishart(const gsl_rng *r, double nu, gsl_matrix *scale,
        gsl_matrix *rwish);
void ran_invwishart(const gsl_rng *r, double nu, gsl_matrix *scale,
        gsl_matrix *riwish);

__END_DECLS

#endif
