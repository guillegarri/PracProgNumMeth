#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <assert.h>


int firststep = 1;

typedef struct PrintInfo {
  int printYN;
  char filename[50];
} pInfo;

void Firststep() {
  firststep = 1;
}

void rkstep32_FSAL(double t,
              double h,
              gsl_vector * yt,
              void f(double t, gsl_vector * y, gsl_vector * dydt),
              gsl_vector * yth,
              gsl_vector * err,
              gsl_vector * FSAL) {
  int n = yt->size;
  int i;

  gsl_vector * k1 = gsl_vector_calloc(n);
  gsl_vector * k2 = gsl_vector_calloc(n);
  gsl_vector * k3 = gsl_vector_calloc(n);
  gsl_vector * ytemp = gsl_vector_calloc(n);

  if (firststep == 1) {
    f(t,yt,k1);
  } else {
    gsl_vector_memcpy(k1, FSAL);
  }


  for (i = 0; i < n; i++) {
    gsl_vector_set(ytemp, i,
                  gsl_vector_get(yt,i) + 0.5*gsl_vector_get(k1,i)*h);
  }

  f(t+0.5*h, ytemp, k2);
  for (i = 0; i < n; i++) {
    gsl_vector_set(ytemp, i,
                  gsl_vector_get(yt,i) + 0.75*gsl_vector_get(k2,i)*h);
  }

  f(t+0.75*h, ytemp, k3);
  for (i = 0; i < n; i++) {
    gsl_vector_set(yth, i,
                  gsl_vector_get(yt,i)
                  + ((2.0/9)*gsl_vector_get(k1,i)
                  + (1.0/3)*gsl_vector_get(k2,i)
                  + (4.0/9)*gsl_vector_get(k3,i))*h);
  }

  f(t+h, yth, FSAL);for (i = 0; i < n; i++) {
    gsl_vector_set(ytemp, i,
                  gsl_vector_get(yt,i)
                  + ((7.0/24)*gsl_vector_get(k1,i)
                  + (1.0/4)*gsl_vector_get(k2,i)
                  + (1.0/3)*gsl_vector_get(k3,i)
                  + (1.0/8)*gsl_vector_get(FSAL,i))*h);

    gsl_vector_set(err,i,gsl_vector_get(yth,i)-gsl_vector_get(ytemp,i));
  }

  firststep = 0;
  gsl_vector_free(k1);
  gsl_vector_free(k2);
  gsl_vector_free(k3);
  gsl_vector_free(ytemp);
}

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
                  pInfo * info
                  ) {
  Firststep();
  int n = yt->size;
  int i;
  double x = *t;
  double a = x;
  double error, ysize, tol;
  double hint = *h; //h internal

  FILE * flptr;
  int printYN = info->printYN;

  //gsl_vector * y = gsl_vector_calloc(n);
  gsl_vector * yh = gsl_vector_calloc(n);
  gsl_vector * dy = gsl_vector_calloc(n);
  gsl_vector * fsal = gsl_vector_calloc(n);

  //gsl_vector_memcpy(y,yt);

  if (printYN>0) {
    flptr = fopen(info->filename,"w");
    fprintf(flptr, "%lg ", x);
    for ( i = 0; i < n; i++) {
      fprintf(flptr, "%lg ", gsl_vector_get(yt,i));
    }
    fprintf(flptr, "\n");
  }

  while (x < b) {
    if (x+hint>b) {hint=b-x;}
    stepper(x,hint,yt,f,yh,dy,fsal);

    error = gsl_blas_dnrm2(dy);
    ysize = gsl_blas_dnrm2(yh);
    tol = (ysize*eps+acc)*sqrt(hint/(b-a));

    if (error<tol) {
      x=x+hint;
      gsl_vector_memcpy(yt,yh);

      if (printYN>0) {
        fprintf(flptr, "%lg ", x);
        for ( i = 0; i < n; i++) {
          fprintf(flptr, "%lg ", gsl_vector_get(yt,i));
        }
        fprintf(flptr, "\n");
      }
    }
    if (error > 0) {
      hint*=pow(tol/error,0.25)*0.95;
    } else {
      hint *= 2;
    }
  }
  *h = hint;

  if (printYN>0) {
    fclose(flptr);
  }


  gsl_vector_free(yh);
  gsl_vector_free(dy);
  gsl_vector_free(fsal);
}

void integral_FSAL(double * x, //start at a
                   double b,
                   double * h,
                   gsl_vector * y,
                   double acc,
                   double eps,
                   void f(double t, gsl_vector * y, gsl_vector * dydt),
                   pInfo * info) {
  assert(y->size==1);
  gsl_vector_set(y,0,0);

  driver_FSAL(x, b, h, y, acc, eps, rkstep32_FSAL, f, info);

}
