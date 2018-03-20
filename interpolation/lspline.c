#include "math.h"
#include "stdio.h"
int search(int n, double * x, double z);

double linterp(int n, double *x, double *y, double z){
  fprintf(stderr, "running linterp\n");
  int i = search(n, x, z);
  double b = (y[i+1] - y[i])/(x[i+1] - x[i]);
  double deltax = z - x[i];
  double value = y[i] + b*deltax;
  return value;
}

double linterp_integ(int n, double *x, double *y, double z){
  int endint = search(n, x, z);
  double integ = 0;
  for (int i = 0; i < endint; i++) {
    double deltax = x[i+1]-x[i];
    double b = (y[i+1] - y[i])/deltax;
    integ += y[i]*deltax + 0.5*b*deltax*deltax;
  }
  double deltax = z - x[endint];
  double b = (y[endint+1] - y[endint])/(x[endint+1] - x[endint]);
  integ += y[endint]*deltax + 0.5*b*deltax*deltax;
  return integ;
}
