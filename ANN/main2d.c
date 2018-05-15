#include<math.h>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

typedef struct {
	int n;
	double (*f)(double,double);
	gsl_vector* data;
} ann2d;

ann2d* ann2d_alloc(int n,double(*f)(double, double));
void ann2d_free(ann2d* nw);
double ann2d_feed(ann2d*nw,double x, double y);
void ann2d_train(ann2d * nw, gsl_vector * X, gsl_vector * Y, gsl_matrix * Z);

double Gwavelet2d(double x, double y) {
  return x*exp(-x*x)*y*exp(-y*y);
}
double functofit(double x, double y) {
  return 6*x + 5*y*y - 8*x*y;
}

int main(int argc, char const *argv[]) {
  double ax=-2, bx=3;
  double ay=-2, by=3;
  int nx = 25;
  int ny = 25;
  int nnodes = 12;
  gsl_vector * X = gsl_vector_alloc(nx);
  gsl_vector * Y = gsl_vector_alloc(ny);
  gsl_matrix * Z = gsl_matrix_alloc(nx,ny);

  for (int i = 0; i < nx; i++) {
	  double xi=ax+(bx-ax)*i/(nx-1);
    gsl_vector_set(X, i, xi);
  }

  for (int j = 0; j < ny; j++) {
    double yj=ay+(by-ay)*j/(ny-1);
    gsl_vector_set(Y, j, yj);
    for (int i = 0; i < nx; i++) {
      double zij = functofit(gsl_vector_get(X,i), yj);
      gsl_matrix_set(Z, i, j, zij);
    }
  }

  ann2d * ANN = ann2d_alloc(nnodes, Gwavelet2d);

  for (int i = 0; i < nnodes; i++) {
    gsl_vector_set(ANN->data, i, ax+(bx-ax)*i/(nnodes-1));
    gsl_vector_set(ANN->data, nnodes+i, 1);
    gsl_vector_set(ANN->data, 2*nnodes+i, ay+(by-ay)*i/(nnodes-1));
    gsl_vector_set(ANN->data, 3*nnodes+i, 1);
    gsl_vector_set(ANN->data, 4*nnodes+i, 1);
  }

  ann2d_train(ANN, X, Y, Z);

  double dx=(bx-ax)/128;
  double dy=(by-ay)/128;
  for (double x=ax;x<=bx;x+=dx) {
    for (double y = ay ; y <= by; y+=dy) {
      double fit = ann2d_feed(ANN,x,y);
      fprintf(stdout, "%g %g %g\n", x, y, fit);
      fprintf(stderr, "%g %g %g\n", x, y, functofit(x,y) );
    }

  }

  ann2d_free(ANN);
  gsl_vector_free(X);
  gsl_vector_free(Y);
  gsl_matrix_free(Z);

  return 0;
}
