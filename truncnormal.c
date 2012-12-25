#include <math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sys.h>

#include "truncnormal.h"

static double ran_left_utruncnormal(const gsl_rng *, double);
static double ran_left_utruncnormal_pos(const gsl_rng *, double);
static double ran_left_utruncnormal_neg(const gsl_rng *, double);

/**
 * Sample from truncated normal.
 */
double
ran_truncnormal(const gsl_rng *r, double a, double b, double m, double s)
{
    double y;
    double am = (a - m) / s;
    double bm = (b - m) / s;

    y = ran_utruncnormal(r, am, bm);

    return s * (y + m);
}

/**
 * Sample from uniform truncated normal.
 */
double
ran_utruncnormal(const gsl_rng *r, double a, double b)
{
    double y;

    if (gsl_isinf(a) && gsl_isinf(b)) {
        y = gsl_ran_ugaussian(r);
    } else if (gsl_isinf(a)) {
        y = -ran_left_utruncnormal(r, -b);
    } else if (gsl_isinf(b)) {
        y = ran_left_utruncnormal(r, a);
    } else {
        do {
            y = gsl_ran_ugaussian(r);
        } while (y < a || y > b);
    }

    return y;
}

static double
ran_left_utruncnormal(const gsl_rng *r, double a)
{
    return (a > 0) ? ran_left_utruncnormal_pos(r, a) :
                     ran_left_utruncnormal_neg(r, a);
}

/* Lower bound is less than or equal to zero.  This is the easy case and we
 * accept-reject naively.
 */
static double
ran_left_utruncnormal_neg(const gsl_rng *r, double a)
{
    double y;

    do {
        y = gsl_ran_ugaussian(r);
    } while (y < a);

    return y;
}

/* See Robert (1992) for details. */
static double
ran_left_utruncnormal_pos(const gsl_rng *r, double a)
{
    double astar;
    double z;
    double u;
    double rho;

    astar = (a + sqrt((a * a) + 4.0)) / 2.0;

    do {

        z = gsl_ran_exponential(r, 1.0 / astar) + a;
        rho = exp(-pow(z - astar, 2) / 2.0);
        u = gsl_ran_flat(r, 0.0, 1.0);
    } while (u > rho);

    return z;
}
