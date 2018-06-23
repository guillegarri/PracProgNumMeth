#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void printM(gsl_matrix * A) {
  for (int i = 0; i < A->size1; i++) {
    for (int j = 0; j < A->size2; j++) {
      fprintf(stderr,"%+-5.3f ",gsl_matrix_get(A,i,j) );
    }
    printf("\n");
  }
}

typedef struct {
	gsl_vector* X;
  gsl_vector* Y;
  gsl_matrix* a;
  gsl_matrix* b;
  gsl_matrix* c;
  gsl_matrix* d;
} bilin;

bilin * bilin_alloc(gsl_vector * X, gsl_vector *Y, gsl_matrix * F);
double bilin2_interp(double x, double y, bilin * sys);
void bilin_free(bilin * sys);

//double rnd(){
//  return ((double)rand()/RAND_MAX);
//}


double dgauss(double x, double y){
  return exp(-x*x)*exp(-y*y);
	//return 1+2*x+3*y+4*x*y;
}


int main(int argc, char const *argv[]) {

  int Nx = 10;
  int Ny = 10;

	int nx = 100;
	int ny = 100;

  double xmin = -3;
  double xmax = 3;
  double ymin = -3;
  double ymax = 3;

  gsl_vector * X = gsl_vector_calloc(Nx);
  gsl_vector * Y = gsl_vector_calloc(Ny);

  gsl_matrix * F = gsl_matrix_calloc(Nx,Ny);

  FILE * truegrid = fopen("true2.out","w");
  FILE * intgrid = fopen("int2.out","w");

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


 bilin * sys = bilin_alloc(X, Y, F);


 for (int i = 0; i < nx; i++) {
    double x = Xmin + i*(Xmax-Xmin)/nx;
    for (int j = 0; j < ny; j++) {
      double y = Ymin + j*(Ymax-Ymin)/ny;

      double val_dgauss = dgauss(x,y);

      double dgauss_int = bilin2_interp(x,y, sys);

      fprintf(truegrid, "%g %g %g\n",x,y, val_dgauss);
      fprintf(intgrid, "%g %g %g\n",x,y, dgauss_int);

    }
    fprintf(truegrid, "\n");
    fprintf(intgrid, "\n");
  }


 bilin_free(sys);

 fclose(truegrid);
 fclose(intgrid);



  gsl_vector_free(X);
  gsl_vector_free(Y);
  gsl_matrix_free(F);

  return 0;
}
