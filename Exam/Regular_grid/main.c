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

double rosenbrock(double x, double y){
  return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}

double himmelblau(double x, double y){
  return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
  //return exp(x)*exp(y);
}

double dgauss(double x, double y){
  return exp(-x*x)*exp(-y*y);
}

double linear(double x, double y){
  return x+y;
}

int main(int argc, char const *argv[]) {
  int Nx = 5;
  int Ny = 5;

  int nx = 50;
  int ny = 30;

  double xmin = -3;
  double xmax = 3;
  double ymin = -3;
  double ymax = 3;

  gsl_vector * X = gsl_vector_calloc(Nx);
  gsl_vector * Y = gsl_vector_calloc(Ny);

  gsl_matrix * Frosenbrock = gsl_matrix_calloc(Nx,Ny);
  gsl_matrix * Fhimmelblau = gsl_matrix_calloc(Nx,Ny);
  gsl_matrix * Fdgauss = gsl_matrix_calloc(Nx,Ny);
  gsl_matrix * Flinear = gsl_matrix_calloc(Nx,Ny);

  FILE * grid_rosenbrock = fopen("grid_rosenbrock.out", "w");
  FILE * grid_himmelblau = fopen("grid_himmelblau.out", "w");
  FILE * grid_dgauss = fopen("grid_dgauss.out", "w");
  FILE * grid_linear = fopen("grid_linear.out", "w");

  FILE * grid_rosenbrock_int = fopen("grid_rosenbrock_int.out", "w");
  FILE * grid_himmelblau_int = fopen("grid_himmelblau_int.out", "w");
  FILE * grid_dgauss_int = fopen("grid_dgauss_int.out", "w");
  FILE * grid_linear_int = fopen("grid_linear_int.out", "w");

  FILE * XYS = fopen("XYvals.out", "w");

  //FILE * XGRID = fopen("XGRID.out", "w");
  //FILE * YGRID = fopen("YGRID.out", "w");

  for (int i = 0; i < Nx; i++) {
    //gsl_vector_set(X,i, xmin+rnd()*(xmax-xmin));
    gsl_vector_set(X,i, xmin+i*(xmax-xmin)/(Nx-1));
  }

  for (int i = 0; i < Ny; i++) {
    //gsl_vector_set(Y,i, ymin+rnd()*(ymax-ymin));
    gsl_vector_set(Y,i, ymin+i*(ymax-ymin)/(Ny-1));
  }

  gsl_sort_vector(X);
  gsl_sort_vector(Y);

  gsl_vector_fprintf(stderr,X,"%g");
  gsl_vector_fprintf(stderr,Y,"%g");

  double Xmin = gsl_vector_get(X,0);
  double Xmax = gsl_vector_get(X,Nx-1);
  double Ymin = gsl_vector_get(Y,0);
  double Ymax = gsl_vector_get(Y,Ny-1);

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      double x = gsl_vector_get(X,i);
      double y = gsl_vector_get(Y,j);
      gsl_matrix_set(Frosenbrock,i,j,rosenbrock(x,y));
      gsl_matrix_set(Fhimmelblau,i,j,himmelblau(x,y));
      gsl_matrix_set(Fdgauss,i,j,dgauss(x,y));
      gsl_matrix_set(Flinear,i,j,linear(x,y));
      fprintf(XYS, "%g %g\n",x, y);
    }
  }

  //gsl_vector * xs = gsl_vector_calloc(nx);
  //gsl_vector * ys = gsl_vector_calloc(ny);

  printf("x \t\ty \t\tRosen \t\tRosen_int \tHimmel \t\tHimmel_int \tdgauss \t\tdgauss_int\n");

  for (int i = 0; i < nx; i++) {
    double x = Xmin + i*(Xmax-Xmin)/nx;
    //gsl_vector_set(xs,i,x);
    for (int j = 0; j < ny; j++) {
      double y = Ymin + j*(Ymax-Ymin)/ny;
      //gsl_vector_set(ys,i,y);

      double val_rosenbrock = rosenbrock(x,y);
      double val_himmelblau = himmelblau(x,y);
      double val_dgauss = dgauss(x,y);
      double val_linear = linear(x,y);

      double rosenbrock_int = interp_bilin(x,y, X, Y, Frosenbrock);
      double himmelblau_int = interp_bilin(x,y, X, Y, Fhimmelblau);
      double dgauss_int = interp_bilin(x,y, X, Y, Fdgauss);
      double linear_int = interp_bilin(x,y,X,Y,Flinear);

      printf("%-10.4g \t%-10.4g \t%-10.4g \t%-10.4g \t%-10.4g \t%-10.4g \t%-10.4g \t%-10.4g\n", x, y, val_rosenbrock, rosenbrock_int, val_himmelblau, himmelblau_int, val_dgauss, dgauss_int);
      fprintf(grid_rosenbrock, "%g %g %g\n",x,y, val_rosenbrock);
      fprintf(grid_rosenbrock_int, "%g %g %g\n",x,y, rosenbrock_int);

      fprintf(grid_himmelblau, "%g %g %g\n",x,y, val_himmelblau);
      fprintf(grid_himmelblau_int, "%g %g %g\n",x,y, himmelblau_int);

      fprintf(grid_dgauss, "%g %g %g\n",x,y, val_dgauss);
      fprintf(grid_dgauss_int, "%g %g %g\n",x,y, dgauss_int);

      fprintf(grid_linear, "%g %g %g\n",x,y, val_linear);
      fprintf(grid_linear_int, "%g %g %g\n",x,y, linear_int);



      //fprintf(XGRID, "%g\n", x);
      //fprintf(YGRID, "%g\n", y);

    }
    fprintf(grid_rosenbrock, "\n");
    fprintf(grid_himmelblau, "\n");
    fprintf(grid_dgauss, "\n");
    fprintf(grid_linear, "\n");

    fprintf(grid_rosenbrock_int, "\n");
    fprintf(grid_himmelblau_int, "\n");
    fprintf(grid_dgauss_int, "\n");
    fprintf(grid_linear_int, "\n");

    //fprintf(XGRID, "\n");
    //fprintf(YGRID, "\n");
  }


  fclose(grid_rosenbrock);
  fclose(grid_himmelblau);
  fclose(grid_dgauss);
  fclose(grid_linear);

  fclose(grid_himmelblau_int);
  fclose(grid_rosenbrock_int);
  fclose(grid_dgauss_int);
  fclose(grid_linear_int);

  fclose(XYS);

  //fclose(XGRID);
  //fclose(YGRID);

  gsl_vector_free(X);
  gsl_vector_free(Y);
  //gsl_vector_free(xs);
  //gsl_vector_free(ys);
  gsl_matrix_free(Frosenbrock);
  gsl_matrix_free(Fhimmelblau);
  gsl_matrix_free(Fdgauss);
  gsl_matrix_free(Flinear);

  return 0;
}
