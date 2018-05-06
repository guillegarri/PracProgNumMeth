#include "math.h"
#include "stdlib.h"

double rnd(){
  return ((double)rand()/RAND_MAX);
}

void randpoint(int dims, double *a, double *b, double *x) {
  for (int i = 0; i < dims; i++) {
    x[i] = a[i] + rnd()* (b[i]-a[i]);
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

  for (int i = 0; i < Numpoints; i++) {
    randpoint(dims, a, b, x);
    funcval = f(x);
    sum +=funcval;
    sumsq +=funcval*funcval;
  }

  double avg = sum/Numpoints;
  double var = sumsq/Numpoints - avg*avg;

  *result = avg*V;
  *error = sqrt(var/Numpoints)*V;
}
