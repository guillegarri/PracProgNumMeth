#include <assert.h>
#include "stdlib.h"
typedef struct {int n; double * x,  *y,  *b,  *c;} qspline;

qspline * qspline_alloc(int n, double * x, double * y){
  qspline * s = (qspline*) malloc(sizeof(qspline));

  s->n =n;
  s->x = (double*) malloc(n*sizeof(double));
  s->y = (double*) malloc(n*sizeof(double));
  s->b = (double*) malloc((n-1)*sizeof(double));
  s->c = (double*) malloc((n-1)*sizeof(double));

  for (int i = 0; i < n; i++) {
    s->x[i]=x[i];
    s->y[i]=y[i];
  }

  int i;
  double dx[n-1], dydx[n-1];
  for (i = 0; i < n-1; i++) {
    dx[i] = x[i+1]-x[i];
    dydx[i] = (y[i+1] - y[i])/dx[i];
  }

  s->c[0] = 0;
  for (i = 0; i < n-2; i++) {
    s->c[i+1] = (dydx[i+1]-dydx[i]-s->c[i]*dx[i])/dx[i+1];
  }
  s->c[n-2] /= 2;
  for (i = n-3; i >= 0; i--) {
    s->c[i] = (dydx[i+1]-dydx[i]-s->c[i+1]*dx[i+1])/dx[i];
  }
  for (i = 0; i < n-1; i++) {
    s->b[i] = dydx[i]-s->c[i]*dx[i];
  }
  return s;
}

double qspline_eval(qspline * s, double z){
  int n = s->n;
  assert(n>1 && z>=s->x[0] && z<=s->x[n-1]);
  int i = 0, j=n-1;
  while (j-i>1) {
    int m=(i+j)/2;
    if (z>s->x[m]) {
      i=m;
    } else {
      j=m;
    }
  }
  return s->y[i] + s->b[i]*(z-s->x[i]) + s->c[i]*(z-s->x[i])*(z-s->x[i]);
}

double qspline_deriv(qspline * s, double z){
  int n = s->n;
  assert(n>1 && z>=s->x[0] && z<=s->x[n-1]);
  int i = 0, j=n-1;
  while (j-i>1) {
    int m=(i+j)/2;
    if (z>s->x[m]) {
      i=m;
    } else {
      j=m;
    }
  }
  return s->b[i] + 2* s->c[i]*(z-s->x[i]);
}

double qspline_int(qspline * s, double z){
  int n = s->n;
  assert(n>1 && z>=s->x[0] && z<=s->x[n-1]);
  int i = 0, j=n-1;
  while (j-i>1) {
    int m=(i+j)/2;
    if (z>s->x[m]) {
      i=m;
    } else {
      j=m;
    }
  }
  double intsum = 0;
  for (int k = 0; k < i; k++) {
    double deltax = s->x[k+1]-s->x[k];
    intsum += s->y[k]*deltax + 0.5*s->b[k]*deltax*deltax + s->c[k]*deltax*deltax*deltax/3.0;
  }
  intsum += s->y[i]*(z-s->x[i]) + 0.5*s->b[i]*(z-s->x[i])*(z-s->x[i]) + s->c[i]*(z-s->x[i])*(z-s->x[i])*(z-s->x[i])/3.0;
  return intsum;
}

void qspline_free(qspline * s) {
  free(s->x);
  free(s->y);
  free(s->b);
  free(s->c);
  free(s);
}
