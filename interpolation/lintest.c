#include "math.h"
#include "stdio.h"

int search(int n, double * x, double z);
double linterp(int n, double *x, double *y, double z);
double linterp_integ(int n, double *x, double *y, double z);

int main(int argc, char const *argv[]) {
int n=20;
double x[n];
double y[n];

double xfine[10*n];
double yfine[10*n];
double yinterp[10*n];
double yinteg[10*n];
double yinteginterp[10*n];

for (int i = 0; i < n; i++) {
  x[i] = 0.5*i;
  y[i] = x[i]*x[i];
}

for (int i = 0; i < 10*n; i++) {
  xfine[i] = 0.05*i;
  yfine[i] = xfine[i]*xfine[i];
  yinteg[i] = xfine[i]*xfine[i]*xfine[i]/3.0;

  if (xfine[i]<x[n-1]) {
    yinterp[i] = linterp(n, x, y, xfine[i]);
    yinteginterp[i] = linterp_integ(n, x, y, xfine[i]);
    printf("%g %g %g %g %g\n",xfine[i], yfine[i], yinterp[i], yinteg[i], yinteginterp[i]);
  }
}



  return 0;
}
