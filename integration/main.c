#include<math.h>
#include<assert.h>
#include<stdio.h>

double integrator24(double f(double), double a, double b, double acc, double eps, double f2, double f3, double *interr,int nrecs);
double integrator(double f(double), double a, double b, double acc, double eps, double *err);

double invsqrt(double x){
  return 1.0/sqrt(x);
}

double lndvsqrt(double x){
  return log(x)/sqrt(x);
}

double testfunction(double x){
  return 4*sqrt(1-(1-x)*(1-x));
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
  printf("Actual error = %g\n", fabs(Qsqrt-2.0/3));


   printf("InvSqrt integrated from 0 to 1:\n");
   double Qinvsqrt = integrator(invsqrt, a, b, acc, eps, &error);
   printf("Calculated = %g\n", Qinvsqrt);
   printf("Exact = %g\n", 2.0);
   printf("Error estimate = %g\n", acc+fabs(Qinvsqrt)*eps);
   printf("Actual error = %g\n", fabs(Qinvsqrt-2.0));

   printf("lndvsqrt integrated from 0 to 1:\n");
   double Qlndvsqrt = integrator(lndvsqrt, a, b, acc, eps, &error);
   printf("Calculated = %g\n", Qlndvsqrt);
   printf("Exact = %g\n", -4.0);
   printf("Error estimate = %g\n", acc+fabs(Qlndvsqrt)*eps);
   printf("Actual error = %g\n", fabs(Qlndvsqrt+4.0));

  printf("testfunction integrated from 0 to 1:\n");
  double Qtest = integrator(testfunction, a, b, acc, eps, &error);
  printf("Calculated = %g\n", Qtest);
  printf("Exact = %g\n", M_PI);
  printf("Error estimate = %g\n", acc+fabs(Qtest)*eps);
  printf("Internal error = %g\n", error);
  printf("Actual error = %g\n", fabs(Qtest-M_PI));



  return 0;
}
