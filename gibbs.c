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
    size_t col;     /* Column (zero-indexed) - day */
};

static double
impmissvar(const gsl_matrix *gsig, const gsl_matrix *gsiginv, size_t col,
        gsl_vector *prod)
{
    gsl_matrix *cov11;
    gsl_vector *cov01;
    double ss;

    cov11 = gsl_matrix_alloc(gsig->size1 - 1, gsig->size2 - 1);
    cov01 = gsl_vector_alloc(gsig->size1 - 1);

    /* prod <- cov01' * inv(cov11) */
    matrix_invert_remove_rowcol(gsiginv, col, cov11, cov01);
    matrix_row_remove_elem(gsig, col, col, cov01);
    gsl_blas_dgemv(CblasTrans, 1.0, cov11, cov01, 0.0, prod);

    /* ss <- cov01' * inv(cov11) * cov01 */
    gsl_blas_ddot(prod, cov01, &ss);

    gsl_matrix_free(cov11); 
    gsl_vector_free(cov01);
    
    return gsl_matrix_get(gsig, col, col) - ss;
}

static double
impmissmean(const gsl_matrix *gdat, const gsl_vector *gthet,
        size_t row, size_t col, const gsl_vector *prod)
{
    double mnmiss;
    size_t i;

    /* mnmiss <- cov01' * inv(cov11) * (y[-i] - gthet[-i]) */
    mnmiss = 0.0;
    for (i = 0; i < col; i++) {
        mnmiss += (gsl_matrix_get(gdat, row, i) - gsl_vector_get(gthet, i)) *
                gsl_vector_get(prod, i);
    }

    for (i = col + 1; i < gdat->size2; i++) {
        mnmiss += (gsl_matrix_get(gdat, row, i) - gsl_vector_get(gthet, i)) *
                gsl_vector_get(prod, i - 1);
    }
    
    return gsl_vector_get(gthet, col) + mnmiss;
}

void
ymissfun(gsl_matrix *gdat, const gsl_vector *gthet, const gsl_matrix *gsig,
        const struct bdl *bdls, size_t nbdls, const gsl_matrix *gsiginv,
        double minmdl, const gsl_rng *rng)
{
    size_t i;
    gsl_matrix *data_copy;
    double *vars;
    gsl_vector **prods;
    size_t ncons;

    data_copy = gsl_matrix_alloc(gdat->size1, gdat->size2);

    gsl_matrix_memcpy(data_copy, gdat);

    ncons = gsig->size2;
    vars = calloc(ncons, sizeof (double));
    prods = calloc(ncons, sizeof (gsl_vector *));
    for (i = 0; i < ncons; i++) {
        prods[i] = gsl_vector_alloc(gsig->size1 - 1);
        vars[i] = impmissvar(gsig, gsiginv, i, prods[i]);
    }

    for (i = 0; i < nbdls; i++) {
        double mean;
        double var;
        double newmiss;

        mean = impmissmean(data_copy, gthet, bdls[i].row, bdls[i].col,
                prods[bdls[i].col]);
        var = vars[bdls[i].col];
        newmiss = ran_truncnormal(rng, minmdl - 10.0, bdls[i].lim, mean,
                sqrt(var));
        gsl_matrix_set(gdat, bdls[i].row, bdls[i].col, newmiss);
    }

    gsl_matrix_free(data_copy);
    free(vars);

    for (i = 0; i < ncons; i++) {
        gsl_vector_free(prods[i]);
    }
    free(prods);
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

void
thetfun(const gsl_matrix *gdat, gsl_vector *gthet, const gsl_matrix *gsig,
        const gsl_matrix *gsiginv, const gsl_rng *rng)
{
    gsl_matrix *chol;
    gsl_vector *means;

    chol = gsl_matrix_alloc(gsig->size1, gsig->size2);
    means = gsl_vector_alloc(gthet->size);

    gsl_matrix_memcpy(chol, gsiginv);
    gsl_matrix_scale(chol, gsig->size1);
    multivariate_mean(gdat, means);
    gsl_blas_dgemv(CblasNoTrans, 1.0, chol, means, 0.0, gthet);

    matrix_update_const_diag(chol, 1.0, 1e-5);
    gsl_linalg_cholesky_decomp(chol);
    gsl_linalg_cholesky_invert(chol);

    gsl_blas_dgemv(CblasNoTrans, 1.0, chol, gthet, 0.0, means);

    ran_multivariate_normal(rng, means, chol, gthet);

    gsl_matrix_free(chol);
    gsl_vector_free(means);
}

void
sigfun(const gsl_matrix *gdat, const gsl_vector *gthet, gsl_matrix *gsig,
        const gsl_rng *rng)
{
    gsl_matrix *swp;
    gsl_matrix *scale;
    double nu;

    swp = gsl_matrix_alloc(gdat->size1, gdat->size2);
    scale = gsl_matrix_alloc(gdat->size2, gdat->size2);

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

    gsl_matrix_free(swp);
    gsl_matrix_free(scale);
}

void
gibbsfun(gsl_matrix *gdat, gsl_vector *gthet, gsl_matrix *gsig,
        const struct bdl *bdls, size_t nbdls, double minmdl, const gsl_rng *rng)
{
    gsl_matrix *gsiginv;

    gsiginv = gsl_matrix_alloc(gsig->size1, gsig->size2);

    gsl_matrix_memcpy(gsiginv, gsig);
    gsl_linalg_cholesky_decomp(gsiginv);
    gsl_linalg_cholesky_invert(gsiginv);

    ymissfun(gdat, gthet, gsig, bdls, nbdls, gsiginv, minmdl, rng);
    thetfun(gdat, gthet, gsig, gsiginv, rng);
    sigfun(gdat, gthet, gsig, rng);

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

void
impute_data(const gsl_matrix *data, const gsl_matrix *mdls,
        const char *output_directory, size_t iterations, size_t burn,
        size_t draws, size_t progress, long seed)
{
    gsl_rng *rng;
    gsl_matrix *gdat;
    gsl_vector *gthet;
    gsl_matrix *gsig;
    gsl_vector *tmp;
    gsl_matrix *imputed_data;
    gsl_vector *mean_mean;
    gsl_matrix *mean_covariance;
    struct bdl *bdls;
    size_t nbdls;
    size_t i;
    size_t j;
    size_t k;
    double minmdl;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, seed);

    gdat = gsl_matrix_alloc(data->size1, data->size2);
    gthet = gsl_vector_alloc(data->size2);
    gsig = gsl_matrix_alloc(gdat->size2, gdat->size2);
    tmp = gsl_vector_alloc(data->size2);

    nbdls = count_bdls(data, mdls);
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

    imputed_data = gsl_matrix_alloc(iterations - burn, nbdls);
    mean_mean = gsl_vector_alloc(data->size2);
    mean_covariance = gsl_matrix_alloc(data->size2, data->size2);

    multivariate_mean(gdat, gthet);
    multivariate_covariance(gdat, gthet, gsig, tmp);
    minmdl = log(gsl_matrix_min(mdls));

    gsl_vector_set_zero(mean_mean);
    gsl_matrix_set_zero(mean_covariance);

    for (i = 0; i < iterations; i++) {
        gibbsfun(gdat, gthet, gsig, bdls, nbdls, minmdl, rng);

        if (i >= burn) {
            gsl_vector_view view = gsl_matrix_row(imputed_data, i - burn);
            size_t j;

            for (j = 0; j < nbdls; j++) {
                gsl_vector_set(&view.vector, j,
                        gsl_matrix_get(gdat, bdls[j].row, bdls[j].col));
            }
            
            gsl_vector_add(mean_mean, gthet);
            gsl_matrix_add(mean_covariance, gsig);
        }

        if (progress > 0 && i % progress == 0) {
            fprintf(stderr, "%ld\n", i);
        }
    }

    gsl_vector_scale(mean_mean, iterations - burn);
    gsl_matrix_scale(mean_covariance, iterations - burn);

    gsl_matrix_free(gdat);
    gsl_vector_free(gthet);
    gsl_matrix_free(gsig);
    gsl_vector_free(tmp);
    gsl_matrix_free(imputed_data);
    gsl_vector_free(mean_mean);
    gsl_matrix_free(mean_covariance);
}
