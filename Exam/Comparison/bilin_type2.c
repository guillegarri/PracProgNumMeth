#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

typedef struct {
	gsl_vector* X;
  gsl_vector* Y;
  gsl_matrix* a;
  gsl_matrix* b;
  gsl_matrix* c;
  gsl_matrix* d;
} bilin;

bilin * bilin_alloc(gsl_vector * X, gsl_vector *Y, gsl_matrix * F){
  int nx = X->size;
  int ny = Y->size;
  assert(nx == F->size1 && ny == F->size2);


  bilin * bilin_sys = malloc(sizeof(bilin));

  gsl_vector * Xcopy = gsl_vector_calloc(nx);
  gsl_vector * Ycopy = gsl_vector_calloc(ny);

  gsl_vector_memcpy(Xcopy, X);
  gsl_vector_memcpy(Ycopy, Y);

  bilin_sys->X=Xcopy;
  bilin_sys->Y=Ycopy;



  gsl_matrix * a = gsl_matrix_calloc(nx-1,ny-1);
  gsl_matrix * b = gsl_matrix_calloc(nx-1,ny-1);
  gsl_matrix * c = gsl_matrix_calloc(nx-1,ny-1);
  gsl_matrix * d = gsl_matrix_calloc(nx-1,ny-1);

  for (int i = 0; i < nx-1; i++) {
    for (int j = 0; j < ny-1; j++) {
      double x1 = gsl_vector_get(X,i);
      double x2 = gsl_vector_get(X,i+1);

      double y1 = gsl_vector_get(Y,j);
      double y2 = gsl_vector_get(Y,j+1);

      double F11 = gsl_matrix_get(F, i,j);
      double F12 = gsl_matrix_get(F, i,j+1);
      double F21 = gsl_matrix_get(F, i+1,j);
      double F22 = gsl_matrix_get(F, i+1,j+1);

      //fprintf(stderr, "%g %g %g %g\n", F11, F12, F21, F22);
      //fprintf(stderr, "%g %g %g %g\n\n", f(x1,y1), f(x1,y2), f(x2,y1), f(x2,y2));

      double aval = (y1*x1*F22-y1*x2*F12-x1*y2*F21+x2*y2*F11)/(x2*y2-x1*y2-x2*y1+x1*y1);
      double bval = -(y1*F22-y1*F12-y2*F21+y2*F11)/((x1-x2)*(y1-y2));
      double cval = -(x1*F22-x1*F21-x2*F12+x2*F11)/((x1-x2)*(y1-y2));
      double dval = (F22 - F21 - F12 + F11)/(x2*y2-x1*y2-x2*y1+x1*y1);

      //fprintf(stderr, "%g %g %g %g\n", aval, bval, cval, dval);

      gsl_matrix_set(a,i,j,aval);
      gsl_matrix_set(b,i,j,bval);
      gsl_matrix_set(c,i,j,cval);
      gsl_matrix_set(d,i,j,dval);
    }
  }

  bilin_sys->a=a;
  bilin_sys->b=b;
  bilin_sys->c=c;
  bilin_sys->d=d;

  return bilin_sys;
}

double bilin2_interp(double x, double y, bilin * sys){
  int nx = sys->X->size;
  int ny = sys->Y->size;
  assert(x >= gsl_vector_get(sys->X,0) && x <= gsl_vector_get(sys->X,nx-1));
  assert(y >= gsl_vector_get(sys->Y,0) && y <= gsl_vector_get(sys->Y,ny-1));

  int i = 0, j=nx-1;
  while (j-i>1) {
    int m=(i+j)/2;
    if (x>gsl_vector_get(sys->X,m)) {
      i=m;
    } else {
      j=m;
    }
  }

  int xi = i;

  i = 0, j=ny-1;
  while (j-i>1) {
    int m=(i+j)/2;
    if (y>gsl_vector_get(sys->Y,m)) {
      i=m;
    } else {
      j=m;
    }
  }

  int yi = i;

  return gsl_matrix_get(sys->a,xi,yi) + x*gsl_matrix_get(sys->b,xi,yi) + y*gsl_matrix_get(sys->c,xi,yi) + x*y*gsl_matrix_get(sys->d,xi,yi);

}

void bilin_free(bilin * sys){
  gsl_vector_free(sys->X);
  gsl_vector_free(sys->Y);
  gsl_matrix_free(sys->a);
  gsl_matrix_free(sys->b);
  gsl_matrix_free(sys->c);
  gsl_matrix_free(sys->d);
  free(sys);
}
