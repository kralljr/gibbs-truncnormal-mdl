#ifndef __TRUNCNORMAL_H__
#define __TRUNCNORMAL_H__

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

double ran_truncnormal(const gsl_rng *r, double a, double b, double m,
        double s);

__END_DECLS

#endif
