#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"
typedef struct {int n; double * x, *y, *b, *c;} cspline;

cspline * cspline_alloc(int n, double * x, double * y);
double cspline_eval(cspline * s, double z);
double cspline_deriv(cspline * s, double z);
double cspline_int(cspline * s, double z);
void cspline_free(cspline * s);

int main(int argc, char const *argv[]) {
  int n=20;
  double x[n];
  double y[n];

  double xfine[10*n];
  double yfine[10*n];
  double yinterp[10*n];
  double yinteg[10*n];
  double yinteginterp[10*n];
  double yderiv[10*n];
  double yderivinterp[10*n];

  for (int i = 0; i < n; i++) {
    x[i] = 0.5*i;
    y[i] = x[i]*x[i];
  }

  cspline * s = cspline_alloc(n, x, y);

  for (int i = 0; i < 10*n; i++) {
    xfine[i] = 0.05*i;
    yfine[i] = xfine[i]*xfine[i];
    yinteg[i] = xfine[i]*xfine[i]*xfine[i]/3.0;
    yderiv[i] = 2*xfine[i];

    if (xfine[i]<x[n-1]) {
      yinterp[i] = cspline_eval(s, xfine[i]);
      yinteginterp[i] = cspline_int(s, xfine[i]);
      yderivinterp[i] = cspline_deriv(s, xfine[i]);
      printf("%g %g %g %g %g %g %g\n",xfine[i], yfine[i], yinterp[i], yinteg[i], yinteginterp[i], yderiv[i], yderivinterp[i]);
    }
  }
  cspline_free(s);

  return 0;
}
