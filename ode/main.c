#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <assert.h>
#include <string.h>

// Useful function and struct
typedef struct PrintInfo {
  int printYN;
  char filename[50];
} pInfo;

pInfo printdef(int YN, char * filename){
  pInfo info;
  info.printYN = YN;
  strcpy(info.filename, filename);
  return info;
}

//Declaring functions from other file

void rkstep32_FSAL(double t,
              double h,
              gsl_vector * yt,
              void f(double t, gsl_vector * y, gsl_vector * dydt),
              gsl_vector * yth,
              gsl_vector * err,
              gsl_vector * FSAL);

void driver_FSAL(
                  double * t,
                  double b,
                  double * h,
                  gsl_vector * yt,
                  double acc,
                  double eps,
                  void stepper(double t,
                                double h,
                                gsl_vector * yt,
                                void f(double t, gsl_vector * y, gsl_vector * dydt),
                                gsl_vector * yth,
                                gsl_vector * err,
                                gsl_vector * FSAL),
                  void f(double t, gsl_vector * y, gsl_vector * dydt),
                  pInfo * info);

void integral_FSAL(double * x, //start at a
                   double b,
                   double * h,
                   gsl_vector * y,
                   double acc,
                   double eps,
                   void f(double t, gsl_vector * y, gsl_vector * dydt),
                   pInfo * info);


//Diffeqs to solve

void sinediffeq(double x, gsl_vector * y, gsl_vector * dydt) {
  assert(y->size == 2);
  gsl_vector_set(dydt,0,gsl_vector_get(y,1));
  gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}

void gaussintdiffeq(double x, gsl_vector * y, gsl_vector * dydt) {
  assert(y->size == 1);
  gsl_vector_set(dydt,0,exp(-x*x));
}

//Main

int main(int argc, char const *argv[]) {

  double x0 = 0;
  double b = 6.28;
  double h = 0.05;
  double acc = 0.0005;
  double eps = 0.01;

  //Solving differential equation for sin(x)
  gsl_vector * y = gsl_vector_calloc(2);
  gsl_vector_set(y,0,0);
  gsl_vector_set(y,1,1);

  pInfo sineprint = printdef(1, "sine.out");

  driver_FSAL(&x0, b, &h, y, acc, eps, rkstep32_FSAL, sinediffeq, &sineprint);

  printf("Solving differential equation for sine from 0 to %g\n", b);
  printf("Final value, y=%g dy/dx=%g\n", gsl_vector_get(y,0),gsl_vector_get(y,1));
  printf("Math values: sin(b)=%g cos(b)=%g\n",sin(b),cos(b));
  printf("Acc = %g eps = %g\n",acc,eps);

  //Solving differential for integral of gaussian from -10 to 10
  gsl_vector * integ = gsl_vector_calloc(1);
  gsl_vector_set(integ,0,0);
  pInfo gaussprint = printdef(0, "NONE");

  double a = -10;
  b = 10;
  h = 0.05;

  printf("Solving differential equation for gaussian integral from %g to %g\n",a, b);


  integral_FSAL(&a, b, &h, integ, acc, eps, gaussintdiffeq, &gaussprint);

  printf("Result, integ=%g\n", gsl_vector_get(integ,0));
  printf("For larger intervals, should converge towards %g\n", sqrt(M_PI));
  //printf("Math values: sin(b)=%g cos(b)=%g\n",sin(b),cos(b));
  //printf("Acc = %g eps = %g\n",acc,eps);




  return 0;
}
