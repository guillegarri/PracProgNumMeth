#include "math.h"
#include "stdio.h"
#include<gsl/gsl_sf_airy.h>

int main() {
  for (double i = -2; i < 2.01; i+=0.01) {
    double x = gsl_sf_airy_Ai(i, GSL_PREC_DOUBLE);
    double y = gsl_sf_airy_Bi(i, GSL_PREC_DOUBLE);
    printf("%lg %lg %lg\n",i, x, y);
  }
  return 0;
}
