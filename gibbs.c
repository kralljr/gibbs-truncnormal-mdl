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
#include "truncnormal.h"
#include "wishart.h"

struct guess {
    gsl_matrix *data;
    gsl_vector *mean;
    gsl_matrix *cvar;
};

/**
 * A single below detection limit observation.
 */
struct bdl {
    double lim;     /* Detection limit */
    size_t row;     /* Row (zero-indexed) - day */
    size_t col;     /* Column (zero-indexed) - day */
};

struct problem {
    const gsl_matrix *data;
    const gsl_matrix *mdls;
    const gsl_rng *rng;
};

/**
 * @chol: Cholesky decomposition of gsig[-wh, -wh]
 */
static double
impmissmn(const gsl_vector *datk, const gsl_vector *gthet, size_t wh,
          const gsl_matrix *chol, const gsl_vector *cov01)
{
    /* XXX: pass, don't allocate */
    gsl_vector *yobs_m_mnobs;
    gsl_vector *soln;
    double mnmiss;
    size_t i;
    size_t j;

    yobs_m_mnobs = gsl_vector_alloc(datk->size - 1);
    soln = gsl_vector_alloc(datk->size - 1);

    j = 0;
    for (i = 0; i < wh; i++, j++) {
        gsl_vector_set(yobs_m_mnobs, j,
                gsl_vector_get(datk, i) - gsl_vector_get(gthet, i));
    }

    for (i = wh + 1; i < datk->size; i++, j++) {
        gsl_vector_set(yobs_m_mnobs, j,
                gsl_vector_get(datk, i) - gsl_vector_get(gthet, i));
    }

    gsl_linalg_cholesky_solve(chol, yobs_m_mnobs, soln);

    gsl_blas_ddot(soln, cov01, &mnmiss);

    gsl_vector_free(yobs_m_mnobs);
    gsl_vector_free(soln);

    return mnmiss + gsl_vector_get(gthet, wh);
}

static double
impmissvar(const gsl_matrix *gsig, size_t wh, const gsl_matrix *chol,
        gsl_vector *cov01)
{
    double ss;

    gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, chol, cov01);
    gsl_blas_ddot(cov01, cov01, &ss);

    return gsl_matrix_get(gsig, wh, wh) - ss;
}

static void
impmisssingle(const gsl_vector *datk, const gsl_vector *gthet,
        const gsl_matrix *gsig, size_t wh, double *mean, double *var)
{
    /* XXX: pass, don't allocate */
    gsl_vector *cov01;
    gsl_matrix *cov11;

    cov01 = gsl_vector_alloc(gsig->size1 - 1);
    cov11 = gsl_matrix_alloc(gsig->size1 - 1, gsig->size2 - 1);

    matrix_row_remove_elem(gsig, wh, wh, cov01);
    matrix_remove_rowcol(gsig, wh, wh, cov11);
    gsl_linalg_cholesky_decomp(cov11);

    *mean = impmissmn(datk, gthet, wh, cov11, cov01);
    *var = impmissvar(gsig, wh, cov11, cov01);

    gsl_vector_free(cov01);
    gsl_matrix_free(cov11);
}

void
ran_multivariate_normal_chol(const gsl_rng *rng, const gsl_vector *mean,
        const gsl_matrix *sigmachol, gsl_vector *z)
{
    size_t i;

    for (i = 0; i < z->size; i++) {
        double r;

        r = gsl_ran_gaussian(rng, 1.0);
        gsl_vector_set(z, i, r);
    }

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, sigmachol, z);
    gsl_vector_add(z, mean);
}

void
ran_multivariate_normal(const gsl_rng *rng, const gsl_vector *mean,
        gsl_matrix *sigma, gsl_vector *z)
{
    gsl_linalg_cholesky_decomp(sigma);
    ran_multivariate_normal_chol(rng, mean, sigma, z);
}

void
ymissfun(gsl_matrix *gdat, const gsl_vector *gthet, const gsl_matrix *gsig,
        const struct bdl *bdls, size_t nbdls, const gsl_rng *rng)
{
    size_t i;

    for (i = 0; i < nbdls; i++) {
        double mean;
        double var;
        double newmiss;

        gsl_vector_const_view row = gsl_matrix_const_row(gdat, bdls[i].row);

        impmisssingle(&row.vector, gthet, gsig, bdls[i].col, &mean, &var);
        newmiss = ran_truncnormal(rng, GSL_NEGINF, bdls[i].lim, mean,
                sqrt(var));

        gsl_matrix_set(gdat, bdls[i].row, bdls[i].col, newmiss);
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

/**
 * New guess for mean: normal.
 *
 * @gdat: missing data (ndays x ncons)
 * @gthet: mean (ncons)
 * @gsig: covariance (ncons x ncons)
 * @rng: random number generator
 *
 * On output, gthet contains the new guess for mean.
 */
void
thetfun(const gsl_matrix *gdat, gsl_vector *gthet, const gsl_matrix *gsig,
        const gsl_rng *rng)
{
    gsl_matrix *chol;
    gsl_vector *means;

    /* XXX: pass, don't allocate */
    chol = gsl_matrix_alloc(gsig->size1, gsig->size2);
    means = gsl_vector_alloc(gthet->size);

    gsl_matrix_memcpy(chol, gsig);
    gsl_linalg_cholesky_decomp(chol);
    gsl_linalg_cholesky_invert(chol);

    gsl_matrix_scale(chol, gsig->size1);
    multivariate_mean(gdat, means);
    gsl_blas_dgemv(CblasNoTrans, 1.0, chol, means, 0.0, gthet);

    matrix_update_const_diag(chol, 1.0, 1e-5);
    gsl_linalg_cholesky_decomp(chol);
    gsl_linalg_cholesky_invert(chol);

    gsl_blas_dgemv(CblasNoTrans, 1.0, chol, gthet, 0.0, means);

    ran_multivariate_normal_chol(rng, means, chol, gthet);

    gsl_matrix_free(chol);
    gsl_vector_free(means);
}

void
sigfun(const gsl_matrix *gdat, const gsl_vector *gthet, gsl_matrix *gsig,
    const gsl_rng *rng)
{
    gsl_matrix *swp;
    gsl_matrix *scale;
    gsl_permutation *permutation;
    double nu;

    /* XXX: pass, don't allocate */
    swp = gsl_matrix_alloc(gdat->size1, gdat->size2);
    scale = gsl_matrix_alloc(gdat->size2, gdat->size2);
    permutation = gsl_permutation_alloc(gdat->size2);

    /* swp <- sweep(gdat, 2, gthet) */
    gsl_matrix_memcpy(swp, gdat);
    matrix_column_sweep(swp, gthet);

    /* scale <- crossprod(swp) + diag(ncol(gdat)) */
    gsl_matrix_set_identity(scale);
    /*gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, swp, 1.0, scale);*/
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, swp, swp, 1.0, scale);

    /* gsig <- riwish(nu, scale) */    
    nu = gdat->size1 + gdat->size2 + 1;
    ran_invwishart(rng, nu, scale, permutation, gsig);

    gsl_matrix_free(swp);
    gsl_matrix_free(scale);
    gsl_permutation_free(permutation);
}

void
gibbsfun(gsl_matrix *gdat, gsl_vector *gthet, gsl_matrix *gsig,
        const struct bdl *bdls, size_t nbdls, const gsl_rng *rng)
{
    ymissfun(gdat, gthet, gsig, bdls, nbdls, rng);
    thetfun(gdat, gthet, gsig, rng);
    sigfun(gdat, gthet, gsig, rng);
}

void
mhwithings(gsl_matrix *gdat, gsl_vector *gthet, gsl_matrix *gsig,
        struct bdl *bdls, size_t nbdls, gsl_rng *rng)
{
    size_t i;

    /* FIXME: pass the constants as arguments */
    for (i = 0; i < 100; i++) {
        gibbsfun(gdat, gthet, gsig, bdls, nbdls, rng);
        fprintf(stderr, "%d\n", i);
    }
}

void
impute_data(const gsl_matrix *data, const gsl_matrix *mdls)
{
    gsl_rng *rng;
    gsl_matrix *gdat;
    gsl_vector *gthet;
    gsl_matrix *gsig;
    gsl_vector *tmp;
    struct bdl *bdls;
    size_t nbdls;
    size_t i;
    size_t j;
    size_t k;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gdat = gsl_matrix_alloc(data->size1, data->size2);
    gthet = gsl_vector_alloc(data->size2);
    gsig = gsl_matrix_alloc(gdat->size2, gdat->size2);
    tmp = gsl_vector_alloc(data->size2);

    nbdls = 0;
    for (i = 0; i < data->size1; i++) {
        for (j = 0; j < data->size2; j++) {
            if (gsl_matrix_get(data, i, j) < gsl_matrix_get(mdls, i, j)) {
                nbdls++;
            }
        }
    }

    bdls = calloc(nbdls, sizeof (struct bdl));
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
        }
    }

    multivariate_mean(gdat, gthet);
    multivariate_covariance(gdat, gthet, gsig, tmp);

    mhwithings(gdat, gthet, gsig, bdls, nbdls, rng);

    gsl_matrix_free(gdat);
    gsl_vector_free(gthet);
    gsl_matrix_free(gsig);
}
