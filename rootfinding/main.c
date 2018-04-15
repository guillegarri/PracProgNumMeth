#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>

int functioncalls = 0;

void qrdecomp(gsl_matrix * A, gsl_matrix * R);
void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x);
void newton(
            void f(gsl_vector * X, gsl_vector * FX),
            gsl_vector * xstart,
            double dx,
            double epsilon);

void newtonWJ(
              void f(gsl_vector * X, gsl_vector * FX, gsl_matrix * JACOBIAN),
                gsl_vector * xstart,
                double epsilon);


void F_one(gsl_vector * x, gsl_vector * fx) {
  double A = 10000;
  double val_one = A*gsl_vector_get(x,0)*gsl_vector_get(x,1)-1.0;
  double val_two = exp(-gsl_vector_get(x,0))+exp(-gsl_vector_get(x,1))-1.0-1.0/A;
  gsl_vector_set(fx, 0, val_one);
  gsl_vector_set(fx, 1, val_two);
  functioncalls++;
}

void GradRosenbrock(gsl_vector * X, gsl_vector * FX) {
  double x = gsl_vector_get(X,0);
  double y = gsl_vector_get(X,1);

  gsl_vector_set(FX,0,2*(200*x*x*x - 200*x*y + x - 1));
  gsl_vector_set(FX,1,200*(y-x*x));
  functioncalls++;
}

void GradHimmelblau(gsl_vector * X, gsl_vector * FX) {
  double x = gsl_vector_get(X,0);
  double y = gsl_vector_get(X,1);

  gsl_vector_set(FX,0, 4*x*(x*x+y-11)+2*x+2*y*y-14);
  gsl_vector_set(FX,1, 4*y*(x+y*y- 7)+2*x*x+2*y-22);
  functioncalls++;
}

void F_oneWJ(gsl_vector * x, gsl_vector * fx, gsl_matrix * J) {
  double A = 10000;
  double val_one = A*gsl_vector_get(x,0)*gsl_vector_get(x,1)-1.0;
  double val_two = exp(-gsl_vector_get(x,0))+exp(-gsl_vector_get(x,1))-1.0-1.0/A;
  gsl_vector_set(fx, 0, val_one);
  gsl_vector_set(fx, 1, val_two);

  gsl_matrix_set(J, 0, 0, A*gsl_vector_get(x,1));
  gsl_matrix_set(J, 0, 1, A*gsl_vector_get(x,0));
  gsl_matrix_set(J, 1, 0, -exp(-gsl_vector_get(x,0)));
  gsl_matrix_set(J, 1, 1, -exp(-gsl_vector_get(x,1)));
  functioncalls++;
}

void GradRosenbrockWJ(gsl_vector * X, gsl_vector * FX, gsl_matrix * J) {
  double x = gsl_vector_get(X,0);
  double y = gsl_vector_get(X,1);

  gsl_vector_set(FX,0,2*(200*x*x*x - 200*x*y + x - 1));
  gsl_vector_set(FX,1,200*(y-x*x));

  gsl_matrix_set(J, 0, 0, 2*(600*x*x - 200*y+1));
  gsl_matrix_set(J, 0, 1, -400*y);
  gsl_matrix_set(J, 1, 0, -400*x);
  gsl_matrix_set(J, 1, 1, 200);
  functioncalls++;
}


int main(int argc, char const *argv[]) {
  gsl_vector * x = gsl_vector_calloc(2);
  gsl_vector * fx = gsl_vector_calloc(2);
  gsl_matrix * J = gsl_matrix_calloc(2,2);


  // First function
  gsl_vector_set(x,0,-2);
  gsl_vector_set(x,1,9);


  newton(F_one, x, 1e-6, 1e-3);

  F_one(x,fx);


  printf("First function: \nx=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("functioncalls = %i \n",functioncalls);
  functioncalls=0;



  //Rosenbrock function
  gsl_vector_set(x,0,0.5);
  gsl_vector_set(x,1,1.5);

  newton(GradRosenbrock, x, 1e-3,1e-5);

  GradRosenbrock(x,fx);

  printf("Rosenbrock function:\nx=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("functioncalls = %i \n",functioncalls);
  functioncalls=0;

  //Himmelblau function
  gsl_vector_set(x,0,2.5);
  gsl_vector_set(x,1,2.5);

  newton(GradHimmelblau, x, 1e-3,1e-5);

  GradHimmelblau(x,fx);

  printf("Himmelblau function:\nx=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("functioncalls = %i \n",functioncalls);
  functioncalls=0;

  //First function with Jacobian
  gsl_vector_set(x,0,-2);
  gsl_vector_set(x,1,9);


  newtonWJ(F_oneWJ, x, 1e-3);

  F_oneWJ(x,fx,J);


  printf("First function: \nx=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("functioncalls = %i \n",functioncalls);
  functioncalls=0;

  //Rosenbrock function with Jacobian
  gsl_vector_set(x,0,0.5);
  gsl_vector_set(x,1,1.5);

  newtonWJ(GradRosenbrockWJ, x,1e-5);

  GradRosenbrockWJ(x,fx,J);

  printf("Rosenbrock function:\nx=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("functioncalls = %i \n",functioncalls);
  functioncalls=0;





  gsl_vector_free(x);
  gsl_vector_free(fx);
  gsl_matrix_free(J);
  return 0;
}
