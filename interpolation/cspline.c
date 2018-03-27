#include <assert.h>
#include "stdlib.h"
#include "stdio.h"
typedef struct {int n; double * x, *y, *b, *c, *d;} cspline;

cspline * cspline_alloc(int n, double * x, double * y){
  cspline * s = (cspline*) malloc(sizeof(cspline));

  s->n =n;
  s->x = (double*) malloc(n*sizeof(double));
  s->y = (double*) malloc(n*sizeof(double));
  s->b = (double*) malloc(n*sizeof(double));
  s->c = (double*) malloc((n-1)*sizeof(double));
  s->d = (double*) malloc((n-1)*sizeof(double));

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

  double D[n], Q[n-1], B[n];
  D[0]=2;
  for (int i = 0; i < n-2; i++) {
    D[i+1] = 2*dx[i]/dx[i+1]+2;
  }
  D[n-1]=2;

  Q[0]=1;
  for (int i = 0; i < n-2; i++) {
    Q[i+1] = dx[i]/dx[i+1];
  }

  for (int i = 0; i < n-2; i++) {
    B[i+1] = 3*(dydx[i]+dydx[i+1]*dx[i]/dx[i+1]);
  }
  B[0] = 3*dydx[0];
  B[n-1] = 3*dydx[n-2];

  for (int i = 1; i < n; i++) {
    D[i] -= Q[i-1]/D[i-1];
    B[i] -= B[i-1]/D[i-1];
  }

  s->b[n-1] = B[n-1]/D[n-1];

  for (int i = n-2; i >=0; i--) {
    s->b[i] = (B[i] - Q[i]*s->b[i+1])/D[i];
  }

  for (int i = 0; i < n-1; i++) {
    s->c[i] = (-2*s->b[i] - s->b[i+1] + 3*dydx[i])/dx[i];
    s->d[i] = (s->b[i] + s->b[i+1] - 2*dydx[i])/dx[i]/dx[i];
  }

  return s;
}

double cspline_eval(cspline * s, double z){
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
  return s->y[i] + s->b[i]*(z-s->x[i]) + s->c[i]*(z-s->x[i])*(z-s->x[i]) + s->d[i]*(z-s->x[i])*(z-s->x[i])*(z-s->x[i]);
}

double cspline_deriv(cspline * s, double z){
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
  return s->b[i] + 2* s->c[i]*(z-s->x[i]) + 3* s->d[i]*(z-s->x[i])*(z-s->x[i]);
}

double cspline_int(cspline * s, double z){
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
    intsum += s->y[k]*deltax + 0.5*s->b[k]*deltax*deltax + s->c[k]*deltax*deltax*deltax/3.0 + s->d[k]*deltax*deltax*deltax*deltax/4.0;
  }
  intsum += s->y[i]*(z-s->x[i]) + 0.5*s->b[i]*(z-s->x[i])*(z-s->x[i]) + s->c[i]*(z-s->x[i])*(z-s->x[i])*(z-s->x[i])/3.0 + s->d[i]*(z-s->x[i])*(z-s->x[i])*(z-s->x[i])*(z-s->x[i])/4.0;
  return intsum;
}

void cspline_free(cspline * s) {
  free(s->x);
  free(s->y);
  free(s->b);
  free(s->c);
  free(s->d);
  free(s);
}
