#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"
typedef struct {int n; double * x, *y, *b, *c;} qspline;

qspline * qspline_alloc(int n, double * x, double * y);
double qspline_eval(qspline * s, double z);
double qspline_deriv(qspline * s, double z);
double qspline_int(qspline * s, double z);
void qspline_free(qspline * s);

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

  qspline * s = qspline_alloc(n, x, y);

  for (int i = 0; i < 10*n; i++) {
    xfine[i] = 0.05*i;
    yfine[i] = xfine[i]*xfine[i];
    yinteg[i] = xfine[i]*xfine[i]*xfine[i]/3.0;
    yderiv[i] = 2*xfine[i];

    if (xfine[i]<x[n-1]) {
      yinterp[i] = qspline_eval(s, xfine[i]);
      yinteginterp[i] = qspline_int(s, xfine[i]);
      yderivinterp[i] = qspline_deriv(s, xfine[i]);
      printf("%g %g %g %g %g %g %g\n",xfine[i], yfine[i], yinterp[i], yinteg[i], yinteginterp[i], yderiv[i], yderivinterp[i]);
    }
  }
  qspline_free(s);

  return 0;
}
