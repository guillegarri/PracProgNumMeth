#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double interp_bilin(double x, double y, gsl_vector * X, gsl_vector * Y, gsl_matrix * F);

double rnd(){
  return ((double)rand()/RAND_MAX);
}


double dgauss(double x, double y){
  return exp(-x*x)*exp(-y*y);
}


int main(int argc, char const *argv[]) {

  int Nints = atoi(argv[1]);

  int Nx = 5;
  int Ny = 5;

  double xmin = -3;
  double xmax = 3;
  double ymin = -3;
  double ymax = 3;

  gsl_vector * X = gsl_vector_calloc(Nx);
  gsl_vector * Y = gsl_vector_calloc(Ny);

  gsl_matrix * F = gsl_matrix_calloc(Nx,Ny);

  for (int i = 0; i < Nx; i++) {
    gsl_vector_set(X,i, xmin+i*(xmax-xmin)/(Nx-1));
  }

  for (int i = 0; i < Ny; i++) {
    gsl_vector_set(Y,i, ymin+i*(ymax-ymin)/(Ny-1));
  }

  gsl_sort_vector(X);
  gsl_sort_vector(Y);

  double Xmin = gsl_vector_get(X,0);
  double Xmax = gsl_vector_get(X,Nx-1);
  double Ymin = gsl_vector_get(Y,0);
  double Ymax = gsl_vector_get(Y,Ny-1);

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      double x = gsl_vector_get(X,i);
      double y = gsl_vector_get(Y,j);
      gsl_matrix_set(F,i,j,dgauss(x,y));
    }
  }

 for (int i = 0; i < Nints; i++) {
   double x = Xmin + rnd()*(Xmax-Xmin);
   double y = Ymin + rnd()*(Ymax-Ymin);
   double Int = interp_bilin(x,y, X, Y, F);
   double val = dgauss(x,y);
   printf("%g %g\n", Int, val);
 }

  gsl_vector_free(X);
  gsl_vector_free(Y);
  gsl_matrix_free(F);

  return 0;
}
