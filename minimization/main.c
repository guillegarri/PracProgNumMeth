#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>

int calls = 0;

int newton(double f(gsl_vector * X, gsl_vector * DF, gsl_matrix * HESSIAN),
            gsl_vector * xstart,
            double eps);

double rosenbrock(gsl_vector * X, gsl_vector * DF, gsl_matrix * H){
  double x = gsl_vector_get(X,0);
  double y = gsl_vector_get(X,1);

  gsl_vector_set(DF,0,2*(200*x*x*x - 200*x*y + x - 1));
  gsl_vector_set(DF,1,200*(y-x*x));

  gsl_matrix_set(H, 0, 0, 2.0*(600.0*x*x - 200.0*y+1.0));
  gsl_matrix_set(H, 0, 1, -400*x);
  gsl_matrix_set(H, 1, 0, -400*x);
  gsl_matrix_set(H, 1, 1, 200);

  calls++;

  return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}

double himmelblau(gsl_vector * X, gsl_vector * DF, gsl_matrix * H){
  double x = gsl_vector_get(X,0);
  double y = gsl_vector_get(X,1);

  gsl_vector_set(DF,0, 4*x*(x*x+y-11)+2*x+2*y*y-14);
  gsl_vector_set(DF,1, 4*y*(x+y*y- 7)+2*x*x+2*y-22);

  gsl_matrix_set(H, 0, 0, 4*(x*x+y-11) + 8*x*x +2);
  gsl_matrix_set(H, 0, 1, 4*x+4*y);
  gsl_matrix_set(H, 1, 0, 4*y+4*x);
  gsl_matrix_set(H, 1, 1, 4*(x+y*y-7)+8*y*y+2);

  calls++;

  return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}


int main(int argc, char const *argv[]) {
  gsl_vector * x = gsl_vector_calloc(2);
  gsl_vector * df = gsl_vector_calloc(2);
  gsl_matrix * H = gsl_matrix_calloc(2,2);
  double y;
  int n;

  // Rosenbrock
  printf("Rosenbrock function:\n");
  gsl_vector_set(x,0,0.6);
  gsl_vector_set(x,1,1.4);
  n = newton(rosenbrock,x,1e-5);
  printf("Minimum at x=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  y = rosenbrock(x,df,H);
  printf("With value %g\n", y);
  printf("and gradient (should be 0) =\n");
  gsl_vector_fprintf(stdout,df,"%g");
  printf("steps taken: %i, functioncalls: %i\n\n", n, calls);

  calls = 0;

  // Himmelblau
  printf("Himmelblau function:\n");
  gsl_vector_set(x,0,2.5);
  gsl_vector_set(x,1,2.5);
  n = newton(himmelblau,x,1e-5);
  printf("Minimum at x=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  y = himmelblau(x,df,H);
  printf("With value %g\n", y);
  printf("and gradient (should be 0) =\n");
  gsl_vector_fprintf(stdout,df,"%g");
  printf("steps taken: %i, functioncalls: %i\n\n", n, calls);

  return 0;
}
