#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

double
ran_truncnormal(const gsl_rng *r, double a, double b, double m, double s)
{
    double y;
    double am = (a - m) / s;
    double bm = (b - m) / s;

    do {
        y = gsl_ran_ugaussian(r);
    } while (y < am || y > bm);

    return s*(y + m);
}
