#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>

int calls = 0;

int broyden(double f(gsl_vector * X, gsl_vector * DF),
            gsl_vector * xstart,
            double eps);

double rosenbrock(gsl_vector * X, gsl_vector * DF){
  double x = gsl_vector_get(X,0);
  double y = gsl_vector_get(X,1);

  gsl_vector_set(DF,0,2*(200*x*x*x - 200*x*y + x - 1));
  gsl_vector_set(DF,1,200*(y-x*x));

  calls++;

  return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}

double himmelblau(gsl_vector * X, gsl_vector * DF){
  double x = gsl_vector_get(X,0);
  double y = gsl_vector_get(X,1);

  gsl_vector_set(DF,0, 4*x*(x*x+y-11)+2*x+2*y*y-14);
  gsl_vector_set(DF,1, 4*y*(x+y*y- 7)+2*x*x+2*y-22);

  calls++;

  return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}

double leastsquares(gsl_vector * X){
  double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
  double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
  double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
  int N = sizeof(t)/sizeof(t[0]);

  double A = gsl_vector_get(X,0);
  double T = gsl_vector_get(X,1);
  double B = gsl_vector_get(X,2);

  double f(double t){
    return A*exp(-t/T)+B;
  }

  double sumsq = 0;
  for (int i = 0; i < N; i++) {
    sumsq += pow( (f(t[i])-y[i])/e[i] ,2);
  }

  return sumsq;
}

double lsqfunc(gsl_vector * X, gsl_vector * DF){
  double fx = leastsquares(X);

  double dx = 1e-6;
  for (int i = 0; i < 3; i++) {
    gsl_vector_set(X,i,gsl_vector_get(X,i)+dx);
    gsl_vector_set(DF,i,(leastsquares(X)-fx)/dx);
    gsl_vector_set(X,i,gsl_vector_get(X,i)-dx);
  }
  calls++;
  return leastsquares(X);
}

int main(int argc, char const *argv[]) {
  gsl_vector * x = gsl_vector_calloc(2);
  gsl_vector * df = gsl_vector_calloc(2);

  gsl_vector * ATB = gsl_vector_calloc(3);
  gsl_vector * dATB = gsl_vector_calloc(3);
  double y;
  int n;

  // Rosenbrock
  printf("Rosenbrock function:\n");
  gsl_vector_set(x,0,0.6);
  gsl_vector_set(x,1,1.4);
  n = broyden(rosenbrock,x,1e-5);
  printf("Minimum at x=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  y = rosenbrock(x,df);
  printf("With value %g\n", y);
  printf("and gradient (should be 0) =\n");
  gsl_vector_fprintf(stdout,df,"%g");
  printf("steps taken: %i, functioncalls: %i\n\n", n, calls);

  calls = 0;

  // Himmelblau
  printf("Himmelblau function:\n");
  gsl_vector_set(x,0,2.5);
  gsl_vector_set(x,1,2.5);
  n = broyden(himmelblau,x,1e-5);
  printf("Minimum at x=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  y = himmelblau(x,df);
  printf("With value %g\n", y);
  printf("and gradient (should be 0) =\n");
  gsl_vector_fprintf(stdout,df,"%g");
  printf("steps taken: %i, functioncalls: %i\n\n", n, calls);

  // Leastsquares
  printf("Least squares fitting:\n");
  gsl_vector_set(ATB,0,1);
  gsl_vector_set(ATB,1,1);
  gsl_vector_set(ATB,2,1);
  n = broyden(lsqfunc,ATB,1e-5);
  printf("Minimum at x=\n");
  gsl_vector_fprintf(stdout,ATB,"%g");
  y = lsqfunc(ATB,dATB);
  printf("With value %g\n", y);
  printf("and gradient (should be 0) =\n");
  gsl_vector_fprintf(stdout,df,"%g");
  printf("steps taken: %i, functioncalls: %i\n\n", n, calls);

  return 0;
}
