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

#pragma omp parallel sections
{

  #pragma omp section
  {
  double error;
  double Qsqrt = integrator(sqrt, a, b, acc, eps, &error);
  printf("Sqrt integrated from 0 to 1:\nCalculated = %g\nExact = %g\nError estimate = %g\nActual error = %g\n\n",
          Qsqrt, 2.0/3, acc+fabs(Qsqrt)*eps,fabs(Qsqrt-2.0/3));}

  #pragma omp section
   {double error;
   double Qinvsqrt = integrator(invsqrt, a, b, acc, eps, &error);
   printf("InvSqrt integrated from 0 to 1:\nCalculated = %g\nExact = %g\nError estimate = %g\nActual error = %g\n\n",
          Qinvsqrt,2.0,acc+fabs(Qinvsqrt)*eps,fabs(Qinvsqrt-2.0));}

   #pragma omp section
   {double error;

   double QCCinvsqrt = ClenshawCurtis(invsqrt, a, b, acc, eps, &error);
   printf("InvSqrt integrated with Clenshaw-Curtis:\nCalculated = %g\nExact = %g\nError estimate = %g\nActual error = %g\n\n",
          QCCinvsqrt, 2.0, acc+fabs(QCCinvsqrt)*eps,fabs(QCCinvsqrt-2.0));}

   #pragma omp section
   {double error;
   double Qlndvsqrt = integrator(lndvsqrt, a, b, acc, eps, &error);
   printf("lndvsqrt integrated from 0 to 1:\nCalculated = %g\nExact = %g\nError estimate = %g\nActual error = %g\n\n",
   Qlndvsqrt, -4.0, acc+fabs(Qlndvsqrt)*eps,fabs(Qlndvsqrt+4.0));}


   #pragma omp section
  {double error;
  double Qtest = integrator(testfunction, a, b, acc, eps, &error);
  printf("testfunction integrated from 0 to 1:\nCalculated = %g\nExact = %g\nError estimate = %g\nInternal error = %g\nActual error = %g\n\n",
          Qtest, M_PI,acc+fabs(Qtest)*eps, error, fabs(Qtest-M_PI));}

  #pragma omp section
  {double error;
  double QCCtest = ClenshawCurtis(testfunction, a, b, acc, eps, &error);
  printf("testfunction integrated with Clenshaw-Curtis:\nCalculated = %g\nExact = %g\nError estimate = %g\nInternal error = %g\nActual error = %g\n\n",
          QCCtest, M_PI, acc+fabs(QCCtest)*eps, error ,fabs(QCCtest-M_PI));}

  #pragma omp section
  {double error;
  double Qgauss = infintegrator(gaussian, -INFINITY, INFINITY, acc, eps, &error);
  printf("Gaussian integrated from -inf to inf:\nCalculated = %g\nExact = %g\nError estimate = %g\nInternal error = %g\nActual error = %g\n\n",
          Qgauss, sqrt(M_PI), acc+fabs(Qgauss)*eps, error, fabs(Qgauss-sqrt(M_PI)));}

  #pragma omp section
  {double error;
  double Qgauss = infintegrator(gaussian, 0, INFINITY, acc, eps, &error);
  printf("Gaussian integrated from 0 to inf:\nCalculated = %g\nExact = %g\nError estimate = %g\nInternal error = %g\nActual error = %g\n\n",
          Qgauss, sqrt(M_PI)/2, acc+fabs(Qgauss)*eps, error, fabs(Qgauss-sqrt(M_PI)/2));}

  #pragma omp section
  {double error;
  double Qgauss = infintegrator(gaussian, -INFINITY, 0, acc, eps, &error);
  printf("Gaussian integrated from -inf to 0:\nCalculated = %g\nExact = %g\nError estimate = %g\nInternal error = %g\nActual error = %g\n\n",
          Qgauss, sqrt(M_PI)/2, acc+fabs(Qgauss)*eps, error, fabs(Qgauss-sqrt(M_PI)/2));}

}



  return 0;
}
