#include<math.h>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>

typedef struct {
	int n;
	double (*f)(double);
	gsl_vector* data;
	} ann;

ann* ann_alloc(int n,double(*f)(double));
void ann_free(ann* nw);
double ann_feed(ann*nw,double x);
void ann_train(ann * nw, gsl_vector * X, gsl_vector * Y);

double Gwavelet(double x) {
  return x*exp(-x*x);
}
double fitfunction(double x) {
//  return x*x - 2*x;
	return exp(x);
}



int main(int argc, char const *argv[]) {
  double a=-2, b=3;
  int numpoints = 30;
  int nnodes = 7;
  gsl_vector * x = gsl_vector_alloc(numpoints);
  gsl_vector * y = gsl_vector_alloc(numpoints);

  for (int i = 0; i < numpoints; i++) {
	  double xi=a+(b-a)*i/(numpoints-1);
	  double yi=fitfunction(xi);
    gsl_vector_set(x, i, xi);
    gsl_vector_set(y, i, yi);
  }

  ann * ANN = ann_alloc(nnodes, Gwavelet);

  for (int i = 0; i < nnodes; i++) {
    gsl_vector_set(ANN->data, i, a+(b-a)*i/(nnodes-1));
    gsl_vector_set(ANN->data, nnodes+i, 1);
    gsl_vector_set(ANN->data, 2*nnodes+i, 1);
  }
  //gsl_vector_fprintf(stderr,ANN->data,"%8.3g");

  ann_train(ANN,x,y);
  //gsl_vector_fprintf(stderr,ANN->data,"%8.3g");

	double dz=(b-a)/128;
  for (double z=a;z<=b;z+=dz) {
    double fit = ann_feed(ANN,z);
    printf("%g %g %g\n",z, fit, fitfunction(z));
  }




  return 0;
}
