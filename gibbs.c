#include <math.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "gibbs.h"
#include "matrix.h"
#include "multivariate.h"
#include "mvnormal.h"
#include "truncnormal.h"
#include "wishart.h"

/**
 * A single below detection limit observation.
 */
struct bdl {
    double lim;     /* Detection limit */
    size_t row;     /* Row (zero-indexed) - day */
    size_t col;     /* Column (zero-indexed) - constituent */
};

struct work {
    size_t time;
    size_t cons;

    gsl_matrix *m_time_cons;
    gsl_matrix *m_cons_cons;
    gsl_vector *v_cons;

    gsl_matrix *m_consm1_consm1;
    gsl_vector *v_consm1;

    double *da_cons;
    gsl_vector **va_consm1_cons;
};

static struct work *
work_new(size_t time, size_t cons)
{
    struct work *w;
    size_t i;

    w = malloc(sizeof (struct work));
    w->time = time;
    w->cons = cons;

    w->m_time_cons = gsl_matrix_alloc(time, cons);
    w->m_cons_cons = gsl_matrix_alloc(cons, cons);
    w->v_cons = gsl_vector_alloc(cons);

    w->m_consm1_consm1 = gsl_matrix_alloc(cons - 1, cons - 1);
    w->v_consm1 = gsl_vector_alloc(cons - 1);

    w->da_cons = calloc(cons, sizeof (double));
    w->va_consm1_cons = calloc(cons, sizeof (struct gsl_vector *));
    for (i = 0; i < cons; i++) {
        w->va_consm1_cons[i] = gsl_vector_alloc(cons - 1);
    }

    return w;
}

static void
work_free(struct work *w)
{
    size_t i;

    gsl_matrix_free(w->m_time_cons);
    gsl_matrix_free(w->m_cons_cons);
    gsl_vector_free(w->v_cons);

    gsl_matrix_free(w->m_consm1_consm1);
    gsl_vector_free(w->v_consm1);

    free(w->da_cons);
    for (i = 0; i < w->cons; i++) {
        gsl_vector_free(w->va_consm1_cons[i]);
    }
    free(w->va_consm1_cons);

    free(w);
}

static void
gsl_matrix_exp(gsl_matrix *m)
{
    size_t i;
    size_t j;

    for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
            gsl_matrix_set(m, i, j, exp(gsl_matrix_get(m, i, j)));
        }
    }
}

/**
 * Subtract v from each row of A.
 */
static void
matrix_column_sweep(gsl_matrix *A, const gsl_vector *v)
{
    size_t i;

    for (i = 0; i < A->size1; i++) {
        gsl_vector_view vA = gsl_matrix_row(A, i);
        gsl_vector_sub(&vA.vector, v);
    }
}

static double
impmissvar(const gsl_matrix *gsig, const gsl_matrix *gsiginv, size_t col,
        gsl_matrix *cov11, gsl_vector *cov01, gsl_vector *prod)
{
    double ss;

    /* prod <- cov01' * inv(cov11) */
    matrix_invert_remove_rowcol(gsiginv, col, cov11, cov01);
    matrix_row_remove_elem(gsig, col, col, cov01);
    gsl_blas_dgemv(CblasTrans, 1.0, cov11, cov01, 0.0, prod);

    /* ss <- cov01' * inv(cov11) * cov01 */
    gsl_blas_ddot(prod, cov01, &ss);
    
    return gsl_matrix_get(gsig, col, col) - ss;
}

static double
impmissmean(const gsl_matrix *gdat, const gsl_vector *gthet,
        size_t row, size_t col, const gsl_vector *prod)
{
    double mnmiss;
    size_t i;
    gsl_vector_const_view vA = gsl_matrix_const_row(gdat, row);

    /* We didn't subtract the guessed means from the row like the R-version
     * does before getting into this function. Do that within each for-loop
     * instead.
     *
     * Directly access vector elements instead of using gsl_vector_get because
     * it's much cheaper - at least for the number of times that we need to do
     * it. */

    /* mnmiss <- cov01' * inv(cov11) * (y[-i] - gthet[-i]) */
    mnmiss = gthet->data[col];
    for (i = 0; i < col; i++) {
        mnmiss += (vA.vector.data[i] - gthet->data[i]) * prod->data[i];
    }

    for (i = col + 1; i < gdat->size2; i++) {
        mnmiss += (vA.vector.data[i] - gthet->data[i]) * prod->data[i - 1];
    }
    
    return mnmiss;
}

void
ymissfun(gsl_matrix *gdat, const gsl_vector *gthet, const gsl_matrix *gsig,
        const struct bdl *bdls, size_t nbdls, const gsl_matrix *gsiginv,
        double minmdl, const gsl_rng *rng, struct work *work)
{
    size_t i;
    double *stds = work->da_cons;
    gsl_vector **prods = work->va_consm1_cons;

    for (i = 0; i < gsig->size2; i++) {
        stds[i] = sqrt(impmissvar(gsig, gsiginv, i, work->m_consm1_consm1,
                work->v_consm1, prods[i]));
    }

    for (i = 0; i < nbdls; i++) {
        double mean;
        double std;
        double newmiss;

        mean = impmissmean(gdat, gthet, bdls[i].row, bdls[i].col,
                prods[bdls[i].col]);
        std = stds[bdls[i].col];
        newmiss = ran_truncnormal(rng, minmdl - 10.0, bdls[i].lim, mean, std);
        gsl_matrix_set(gdat, bdls[i].row, bdls[i].col, newmiss);
    }
}

void
thetfun(const gsl_matrix *gdat, gsl_vector *gthet, const gsl_matrix *gsig,
        const gsl_matrix *gsiginv, const gsl_rng *rng, struct work *work)
{
    gsl_matrix *chol = work->m_cons_cons;
    gsl_vector *means = work->v_cons;

    gsl_matrix_memcpy(chol, gsiginv);
    gsl_matrix_scale(chol, gsig->size1);
    multivariate_mean(gdat, means);
    gsl_blas_dgemv(CblasNoTrans, 1.0, chol, means, 0.0, gthet);

    matrix_update_const_diag(chol, 1.0, 1e-5);
    gsl_linalg_cholesky_decomp(chol);
    gsl_linalg_cholesky_invert(chol);

    gsl_blas_dgemv(CblasNoTrans, 1.0, chol, gthet, 0.0, means);

    ran_multivariate_normal(rng, means, chol, gthet);
}

void
sigfun(const gsl_matrix *gdat, const gsl_vector *gthet, gsl_matrix *gsig,
        const gsl_rng *rng, struct work *work)
{
    gsl_matrix *swp = work->m_time_cons;
    gsl_matrix *scale = work->m_cons_cons;
    double nu;

    /* swp <- sweep(gdat, 2, gthet) */
    gsl_matrix_memcpy(swp, gdat);
    matrix_column_sweep(swp, gthet);

    /* scale <- crossprod(swp) + diag(ncol(gdat)) */
    gsl_matrix_set_identity(scale);
    /*gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, swp, 1.0, scale);*/
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, swp, swp, 1.0, scale);

    /* gsig <- riwish(nu, scale) */    
    nu = gdat->size1 + gdat->size2 + 1;
    ran_invwishart(rng, nu, scale, gsig);
}

void
gibbsfun(gsl_matrix *gdat, gsl_vector *gthet, gsl_matrix *gsig,
        const struct bdl *bdls, size_t nbdls, double minmdl, const gsl_rng *rng,
        struct work *work)
{
    gsl_matrix *gsiginv;

    gsiginv = gsl_matrix_alloc(gsig->size1, gsig->size2);

    gsl_matrix_memcpy(gsiginv, gsig);
    gsl_linalg_cholesky_decomp(gsiginv);
    gsl_linalg_cholesky_invert(gsiginv);

    ymissfun(gdat, gthet, gsig, bdls, nbdls, gsiginv, minmdl, rng, work);
    thetfun(gdat, gthet, gsig, gsiginv, rng, work);
    sigfun(gdat, gthet, gsig, rng, work);

    gsl_matrix_free(gsiginv);
}

static size_t
count_bdls(const gsl_matrix *data, const gsl_matrix *mdls)
{
    size_t nbdls;
    size_t i;
    size_t j;

    nbdls = 0;
    for (i = 0; i < data->size1; i++) {
        for (j = 0; j < data->size2; j++) {
            if (gsl_matrix_get(data, i, j) < gsl_matrix_get(mdls, i, j)) {
                nbdls++;
            }
        }
    }
   
    return nbdls; 
}

static double
find_bdls(const gsl_matrix *data, const gsl_matrix *mdls, gsl_matrix *gdat,
        struct bdl *bdls)
{
    double minmdl;
    size_t i;
    size_t j;
    size_t k;

    minmdl = gsl_matrix_get(mdls, 0, 0);
    k = 0;
    for (i = 0; i < data->size1; i++) {
        for (j = 0; j < data->size2; j++) {
            double delem;
            double melem;

            delem = gsl_matrix_get(data, i, j);
            melem = gsl_matrix_get(mdls, i, j); 
            if (delem < melem) {
                gsl_matrix_set(gdat, i, j, log(melem) / 2.0);
                bdls[k].lim = log(melem) / 2.0;
                bdls[k].row = i;
                bdls[k].col = j;
                k++;
            } else {
                gsl_matrix_set(gdat, i, j, log(delem));
            }

            if (melem < minmdl) {
                minmdl = melem;
            }
        }
    }

    return log(minmdl);
}

static int
cmp_size_t(const void *s1, const void *s2)
{
    const size_t x1 = *(const size_t *) s1;
    const size_t x2 = *(const size_t *) s2; 

    if (x1 < x2) {
        return -1;
    } else if (x1 > x2) {
        return 1;
    } else {
        return 0;
    }
}

static void
qsort_size_t(size_t *s, size_t n)
{
    qsort(s, n, sizeof (size_t), cmp_size_t);
}

static void
generate_draws(gsl_rng *rng, size_t *iterations, size_t ndraws, size_t min,
        size_t max)
{
    size_t i = 0;

    while (i < ndraws) {
        size_t j;
        size_t new_draw;
        int repeat;

        new_draw = gsl_rng_uniform_int(rng, max - min) + min;
        repeat = 0;
        for (j = 0; j < i; j++) {
            if (new_draw == iterations[j]) {
                repeat = 1;
                break;
            }
        }

        if (!repeat) {
            iterations[i] = new_draw;
            i++;
        }
    }

    qsort_size_t(iterations, ndraws);
}

static void
gibbs_problem_alloc(struct gibbs_problem *p)
{
    size_t time;
    size_t cons;
    size_t i;

    time = p->data->size1;
    cons = p->data->size2;

    p->dindexes = calloc(p->draws, sizeof (size_t));
    p->ddata = calloc(p->draws, sizeof (gsl_matrix *));
    for (i = 0; i < p->draws; i++) {
        p->ddata[i] = gsl_matrix_alloc(time, cons);
    }
    p->mthet = gsl_vector_alloc(cons);
    p->msig = gsl_matrix_alloc(cons, cons);
}

void
gibbs_problem_exec(struct gibbs_problem *p)
{
    size_t time;
    size_t cons;

    gsl_rng *rng;

    gsl_matrix *gdat;
    gsl_vector *gthet;
    gsl_matrix *gsig;

    size_t nbdls;
    struct bdl *bdls;
    double minmdl;
    struct work *work;

    size_t i;
    size_t j;

    time = p->data->size1;
    cons = p->data->size2;

    gibbs_problem_alloc(p);
    rng = gsl_rng_alloc(gsl_rng_taus);
    work = work_new(time, cons);

    gdat = gsl_matrix_alloc(time, cons);
    gthet = gsl_vector_alloc(cons);
    gsig = gsl_matrix_alloc(cons, cons);

    nbdls = count_bdls(p->data, p->mdls);
    bdls = calloc(nbdls, sizeof (struct bdl));

    minmdl = find_bdls(p->data, p->mdls, gdat, bdls);
    multivariate_mean(gdat, gthet);
    multivariate_covariance(gdat, gthet, gsig, p->mthet);

    gsl_vector_set_zero(p->mthet);
    gsl_matrix_set_zero(p->msig);

    gsl_rng_set(rng, p->seed);
    generate_draws(rng, p->dindexes, p->draws, p->burn, p->iterations);

    j = 0;
    for (i = 0; i < p->iterations; i++) {
        gibbsfun(gdat, gthet, gsig, bdls, nbdls, minmdl, rng, work);

        if (i >= p->burn) {
            gsl_vector_add(p->mthet, gthet);
            gsl_matrix_add(p->msig, gsig);
        }

        if (j < p->draws && i == p->dindexes[j]) {
            gsl_matrix_memcpy(p->ddata[j], gdat);
            gsl_matrix_exp(p->ddata[j]);
            
            j++;
        }
    }

    gsl_vector_scale(p->mthet, 1.0 / (p->iterations - p->burn));
    gsl_matrix_scale(p->msig, 1.0 / (p->iterations - p->burn));

    work_free(work);
    gsl_rng_free(rng);
    gsl_matrix_free(gdat);
    gsl_vector_free(gthet);
    gsl_matrix_free(gsig);
    free(bdls);
}

void
gibbs_problem_free(struct gibbs_problem *p)
{
    size_t i;

    free(p->dindexes);
    for (i = 0; i < p->draws; i++) {
        gsl_matrix_free(p->ddata[i]);
    }
    free(p->ddata);
    gsl_vector_free(p->mthet);
    gsl_matrix_free(p->msig);
}
