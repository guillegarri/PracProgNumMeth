#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double interp_bilin(double x, double y, gsl_vector * X, gsl_vector * Y, gsl_matrix * F);

int main(int argc, char const *argv[]) {
  int Nx = 5;
  int Ny = 5;

  int nx = 50;
  int ny = 50;

  double xmin = -3;
  double xmax = 3;
  double ymin = -3;
  double ymax = 3;

  gsl_vector * X = gsl_vector_calloc(Nx);
  gsl_vector * Y = gsl_vector_calloc(Ny);

  gsl_matrix * F = gsl_matrix_calloc(Nx,Ny);

  

  return 0;
}
