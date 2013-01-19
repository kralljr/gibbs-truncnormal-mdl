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

struct gibbs_problem {
    const gsl_matrix *data;
    const gsl_matrix *mdls;

    size_t iterations;
    size_t burn;
    size_t draws;
    long seed;

    gsl_matrix *dat;
    gsl_vector *mthet;
    gsl_matrix *msig;
    gsl_matrix **ddata;
    size_t *dindexes;
};

void gibbs_problem_exec(struct gibbs_problem *p);
void gibbs_problem_free(struct gibbs_problem *p);


__END_DECLS

#endif
