#include "math.h"
#include "stdlib.h"
#include <gsl/gsl_rng.h>

void randpoint(int dims, double *a, double *b, double *x, gsl_rng *RNG) {
  for (int i = 0; i < dims; i++) {
    x[i] = a[i] + gsl_rng_uniform(RNG)* (b[i]-a[i]);
  }
}

void regularmcMP(int dims, double *a, double *b, double f(double *x), int Numpoints, double * result, double * error) {
  double V=1;
  for (int i = 0; i < dims; i++) {
    V*=b[i]-a[i];
  }

  double sum1 = 0;
  double sumsq1 = 0;
  double funcval1;
  double x1[dims];

  double sum2 = 0;
  double sumsq2 = 0;
  double funcval2;
  double x2[dims];

  gsl_rng * RNG1 = gsl_rng_alloc( gsl_rng_ranlux389);
  gsl_rng * RNG2 = gsl_rng_alloc( gsl_rng_ranlux389);

  gsl_rng_set(RNG1,1);
  gsl_rng_set(RNG2,2);

#pragma omp parallel sections
{

#pragma omp section
{
  for (int i = 0; i < Numpoints/2; i++) {
    randpoint(dims, a, b, x1, RNG1);
    funcval1 = f(x1);
    sum1 +=funcval1;
    sumsq1 +=funcval1*funcval1;
  }
}

#pragma omp section
{
  for (int i = Numpoints/2; i < Numpoints; i++) {
    randpoint(dims, a, b, x2, RNG2);
    funcval2 = f(x2);
    sum2 +=funcval2;
    sumsq2 +=funcval2*funcval2;
  }
}
}

  double avg = (sum1+sum2)/Numpoints;
  double var = (sumsq1+sumsq2)/Numpoints - avg*avg;

  *result = avg*V;
  *error = sqrt(var/Numpoints)*V;

  gsl_rng_free(RNG1);
  gsl_rng_free(RNG2);
}
