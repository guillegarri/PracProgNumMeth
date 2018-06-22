#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

double interp_bilin(double x, double y, gsl_vector * X, gsl_vector * Y, gsl_matrix * F){
  int nx = X->size;
  int ny = Y->size;
  //fprintf(stderr, "nx = %i, ny = %i\n", nx, ny);
  //fprintf(stderr, "xmin = %g x=%g xmax=%g\n", gsl_vector_get(X,0), x, gsl_vector_get(X,nx-1));

  assert(nx == F->size1 && ny == F->size2);
  assert(x >= gsl_vector_get(X,0) && x <= gsl_vector_get(X,nx-1));
  assert(y >= gsl_vector_get(Y,0) && y <= gsl_vector_get(Y,ny-1));

  int i = 0, j=nx-1;
  while (j-i>1) {
    int m=(i+j)/2;
    if (x>gsl_vector_get(X,m)) {
      i=m;
    } else {
      j=m;
    }
  }

  int xi = i;
  double x1 = gsl_vector_get(X,xi);
  double x2 = gsl_vector_get(X,xi+1);

  i = 0, j=ny-1;
  while (j-i>1) {
    int m=(i+j)/2;
    if (y>gsl_vector_get(Y,m)) {
      i=m;
    } else {
      j=m;
    }
  }

  int yi = i;
  double y1 = gsl_vector_get(Y,yi);
  double y2 = gsl_vector_get(Y,yi+1);

  fprintf(stderr, "x1=%g x=%g x2=%g\n", x1, x, x2);
  fprintf(stderr, "y1=%g y=%g y2=%g\n\n", y1, y, y2);



  double F11 = gsl_matrix_get(F, xi,yi);
  double F12 = gsl_matrix_get(F, xi,yi+1);
  double F21 = gsl_matrix_get(F, xi+1,yi);
  double F22 = gsl_matrix_get(F, xi+1,yi+1);

  double a = (y1*x1*F22-y1*x2*F12-x1*y2*F21+x2*y2*F11)/(x2*y2-x1*y2-x2*y1+x1*y1);
  double b = -(y1*F22-y1*F12-y2*F21+y2*F11)/((x1-x2)*(y1-y2));
  double c = -(x1*F22-x1*F21-x2*F12+x2*F11)/((x1-x2)*(y1-y2));
  double d = (F22 - F21 - F12 + F11)/(x2*y2-x1*y2-x2*y1+x1*y1);

  return a + b*x + c*y + d*x*y;
}
