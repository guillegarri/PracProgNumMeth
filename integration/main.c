#include<math.h>
#include<assert.h>
#include<stdio.h>

double integrator24(double f(double), double a, double b, double acc, double eps, double f2, double f3, double *interr,int nrecs);
double integrator(double f(double), double a, double b, double acc, double eps, double *err);
double infintegrator(double f(double), double a, double b, double acc, double eps, double *err);
double ClenshawCurtis(double f(double), double a, double b, double acc, double eps, double *err);

double invsqrt(double x){
  return 1.0/sqrt(x);
}

double lndvsqrt(double x){
  return log(x)/sqrt(x);
}

double testfunction(double x){
  return 4*sqrt(1-(1-x)*(1-x));
}

double gaussian(double x){
  return exp(-x*x);
}

int main(int argc, char const *argv[]) {

  double a = 0;
  double b = 1;
  double acc =1e-8;
  double eps =1e-8;
  double error;
  printf("Sqrt integrated from 0 to 1:\n");
  double Qsqrt = integrator(sqrt, a, b, acc, eps, &error);
  printf("Calculated = %g\n", Qsqrt);
  printf("Exact = %g\n", 2.0/3);
  printf("Error estimate = %g\n", acc+fabs(Qsqrt)*eps);
  printf("Actual error = %g\n\n", fabs(Qsqrt-2.0/3));


   printf("InvSqrt integrated from 0 to 1:\n");
   double Qinvsqrt = integrator(invsqrt, a, b, acc, eps, &error);
   printf("Calculated = %g\n", Qinvsqrt);
   printf("Exact = %g\n", 2.0);
   printf("Error estimate = %g\n", acc+fabs(Qinvsqrt)*eps);
   printf("Actual error = %g\n\n", fabs(Qinvsqrt-2.0));


   printf("InvSqrt integrated with Clenshaw-Curtis:\n");
   double QCCinvsqrt = ClenshawCurtis(invsqrt, a, b, acc, eps, &error);
   printf("Calculated = %g\n", QCCinvsqrt);
   printf("Exact = %g\n", 2.0);
   printf("Error estimate = %g\n", acc+fabs(QCCinvsqrt)*eps);
   printf("Actual error = %g\n\n", fabs(QCCinvsqrt-2.0));

   printf("lndvsqrt integrated from 0 to 1:\n");
   double Qlndvsqrt = integrator(lndvsqrt, a, b, acc, eps, &error);
   printf("Calculated = %g\n", Qlndvsqrt);
   printf("Exact = %g\n", -4.0);
   printf("Error estimate = %g\n", acc+fabs(Qlndvsqrt)*eps);
   printf("Actual error = %g\n\n", fabs(Qlndvsqrt+4.0));



  printf("testfunction integrated from 0 to 1:\n");
  double Qtest = integrator(testfunction, a, b, acc, eps, &error);
  printf("Calculated = %g\n", Qtest);
  printf("Exact = %g\n", M_PI);
  printf("Error estimate = %g\n", acc+fabs(Qtest)*eps);
  printf("Internal error = %g\n", error);
  printf("Actual error = %g\n\n", fabs(Qtest-M_PI));

  printf("testfunction integrated with Clenshaw-Curtis:\n");
  double QCCtest = ClenshawCurtis(testfunction, a, b, acc, eps, &error);
  printf("Calculated = %g\n", QCCtest);
  printf("Exact = %g\n", M_PI);
  printf("Error estimate = %g\n", acc+fabs(QCCtest)*eps);
  printf("Internal error = %g\n", error);
  printf("Actual error = %g\n\n", fabs(QCCtest-M_PI));

  printf("Gaussian integrated from -inf to inf:\n");
  double Qgauss = infintegrator(gaussian, -INFINITY, INFINITY, acc, eps, &error);
  printf("Calculated = %g\n", Qgauss);
  printf("Exact = %g\n", sqrt(M_PI));
  printf("Error estimate = %g\n", acc+fabs(Qgauss)*eps);
  printf("Internal error = %g\n", error);
  printf("Actual error = %g\n\n", fabs(Qgauss-sqrt(M_PI)));

  printf("Gaussian integrated from 0 to inf:\n");
  Qgauss = infintegrator(gaussian, 0, INFINITY, acc, eps, &error);
  printf("Calculated = %g\n", Qgauss);
  printf("Exact = %g\n", sqrt(M_PI)/2);
  printf("Error estimate = %g\n", acc+fabs(Qgauss)*eps);
  printf("Internal error = %g\n", error);
  printf("Actual error = %g\n\n", fabs(Qgauss-sqrt(M_PI)/2));

  printf("Gaussian integrated from -inf to 0:\n");
  Qgauss = infintegrator(gaussian, -INFINITY, 0, acc, eps, &error);
  printf("Calculated = %g\n", Qgauss);
  printf("Exact = %g\n", sqrt(M_PI)/2);
  printf("Error estimate = %g\n", acc+fabs(Qgauss)*eps);
  printf("Internal error = %g\n", error);
  printf("Actual error = %g\n\n", fabs(Qgauss-sqrt(M_PI)/2));





  return 0;
}
