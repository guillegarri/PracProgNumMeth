#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void regularmcMP(int dims, double *a, double *b, double f(double *x), int Numpoints, double * result, double * error);

double twodgauss(double *x){
  return exp(-x[0]*x[0]-x[1]*x[1]);
}

double distsquared3d(double *x){
  return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

double testfunction(double * x){
  return 1.0/((1-cos(x[0])*cos(x[1])*cos(x[2]))*M_PI*M_PI*M_PI);
}

int main(int argc, char const *argv[]) {

  //First test function
  int D = 2;
  int N = 10000;
  double a1[2] = {-10.0, -10.0};
  double b1[2] = {10.0, 10.0};
  double result;
  double error;

  regularmcMP(D, a1, b1, twodgauss, N, &result, &error);
  printf("First test function, a 2d Gaussian\n");
  printf("Result: %g\n", result);
  printf("Target result: %g\n", M_PI);
  printf("Estimated error: %g\n", error);
  printf("Real error: %g\n\n", fabs(result-M_PI));


  //Second test function
  D=3;
  double a2[3] ={-5.0, -5.0, -5.0};
  double b2[3] ={5.0, 5.0, 5.0};

  regularmcMP(D, a2, b2, distsquared3d, N, &result, &error);
  printf("Second test function, 3d distance squared\n");
  printf("Result: %g\n", result);
  printf("Target result: %g\n", 25000.0);
  printf("Estimated error: %g\n", error);
  printf("Real error: %g\n\n", fabs(result-25000.0));

  //Provided test function
  //N*=1000;
  double a3[3] ={0,0,0};
  double b3[3] ={M_PI, M_PI, M_PI};
  regularmcMP(D, a3, b3, testfunction, N, &result, &error);
  printf("Third test function, provided function\n");
  printf("Result: %g\n", result);
  printf("Target result: %g\n", 1.393203929685);
  printf("Estimated error: %g\n", error);
  printf("Real error: %g\n\n", fabs(result-1.393203929685));

  //Error as function of N, using 2d Gaussian
  D=2;
  for (N = 10; N < 1e9; N*=10) {
    regularmcMP(D, a1, b1, twodgauss, N, &result, &error);
    fprintf(stderr, "%i %g %g %g %g\n",N,result, M_PI, error, fabs(result-M_PI) );
  }

  return 0;}
