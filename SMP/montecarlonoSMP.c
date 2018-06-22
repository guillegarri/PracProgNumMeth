#include "math.h"
#include "stdlib.h"
#include <gsl/gsl_rng.h>



void randpoint(int dims, double *a, double *b, double *x, gsl_rng *RNG) {
  for (int i = 0; i < dims; i++) {
    x[i] = a[i] + gsl_rng_uniform(RNG)* (b[i]-a[i]);
  }
}

void regularmc(int dims, double *a, double *b, double f(double *x), int Numpoints, double * result, double * error) {
  double V=1;
  for (int i = 0; i < dims; i++) {
    V*=b[i]-a[i];
  }

  double sum = 0;
  double sumsq = 0;
  double funcval;
  double x[dims];

  gsl_rng * RNG = gsl_rng_alloc( gsl_rng_ranlux389);

  gsl_rng_set(RNG,3);

  for (int i = 0; i < Numpoints; i++) {
    randpoint(dims, a, b, x, RNG);
    funcval = f(x);
    sum +=funcval;
    sumsq +=funcval*funcval;
  }

  double avg = sum/Numpoints;
  double var = sumsq/Numpoints - avg*avg;

  *result = avg*V;
  *error = sqrt(var/Numpoints)*V;
}
