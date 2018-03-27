#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"
typedef struct {int n; double * x, *y, *b, *c;} qspline;
typedef struct {int n; double * x, *y, *b, *c, *d;} cspline;

int search(int n, double * x, double z);
double linterp(int n, double *x, double *y, double z);
double linterp_integ(int n, double *x, double *y, double z);
qspline * qspline_alloc(int n, double * x, double * y);
double qspline_eval(qspline * s, double z);
double qspline_deriv(qspline * s, double z);
double qspline_int(qspline * s, double z);
void qspline_free(qspline * s);
cspline * cspline_alloc(int n, double * x, double * y);
double cspline_eval(cspline * s, double z);
double cspline_deriv(cspline * s, double z);
double cspline_int(cspline * s, double z);
void cspline_free(cspline * s);

int main(int argc, char const *argv[]) {

  int N=6;
  double x[] = {0, 2, 4, 6, 8, 10};
  double y[] = {0, 4, 16, 36, 64, 100};

  qspline * qs = qspline_alloc(N, x, y);
  cspline * cs = cspline_alloc(N, x, y);

  int n=100;



  for (int i = 0; i < n; i++) {
    double z = i/10.0;
   // z ly[n], qy[n], cy[n], lint[n], qint[n], cint[n], qd[n], cd[n];
    printf("%g %g %g %g %g %g %g %g %g\n",
            z,
            linterp(N, x, y, z), qspline_eval(qs, z), cspline_eval(cs, z),
            linterp_integ(N, x, y, z), qspline_int(qs, z), cspline_int(cs, z),
            qspline_deriv(qs, z), cspline_deriv(cs, z));
  }

  qspline_free(qs);
  cspline_free(cs);

  return 0;
}
