#ifndef __GIBBS_H__
#define __GIBBS_H__

#include <stddef.h>

#include <gsl/gsl_matrix.h>

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

void impute_data(const gsl_matrix *data, const gsl_matrix *mdls,
        const char *output_directory, size_t iterations, size_t skip,
        size_t draws, size_t progress, long seed);

__END_DECLS

#endif
