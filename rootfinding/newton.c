#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
//#include <stdio.h>

void qrdecomp(gsl_matrix * A, gsl_matrix * R);
void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x);

void newton(
            void f(gsl_vector * X, gsl_vector * FX),
            gsl_vector * xstart,
            double dx,
            double epsilon) {
  int n = xstart->size;
  int maxsteps = 1000000;
  int stepsdone = 0;
  gsl_matrix * J = gsl_matrix_calloc(n,n);
  gsl_matrix * R = gsl_matrix_calloc(n,n);
  gsl_vector * x = gsl_vector_calloc(n);
  gsl_vector * fx = gsl_vector_calloc(n);
  gsl_vector * mfx = gsl_vector_calloc(n);
  gsl_vector * xplus = gsl_vector_calloc(n);
  gsl_vector * fxplus = gsl_vector_calloc(n);
  gsl_vector * deltax = gsl_vector_calloc(n);
  double lambda;
  double fxsize;
  double fxplussize;

  gsl_vector_memcpy(x, xstart);
  f(x,fx);


  do {
    for (int k = 0; k < n; k++) {
      gsl_vector_memcpy(xplus,x);
      gsl_vector_set(xplus,k,gsl_vector_get(xplus,k)+dx);
      f(xplus, fxplus);
      gsl_vector_sub(fxplus, fx);
      for (int i = 0; i < n; i++) {
        //dfdx = (fxplus - fx)/dx;
        gsl_matrix_set(J,i,k, gsl_vector_get(fxplus,i)/dx);
      }
    }

    gsl_vector_memcpy(mfx, fx);
    gsl_vector_scale(mfx,-1);

    qrdecomp(J,R);
    qrbacksub(J,R,mfx,deltax);

    lambda = 1.0;
    do {
      fxsize = gsl_blas_dnrm2(fx);

      gsl_vector_memcpy(xplus,x);
      gsl_vector_scale(deltax,lambda);
      gsl_vector_add(xplus,deltax);
      gsl_vector_scale(deltax, 1.0/lambda);
      f(xplus,fxplus);

      fxplussize = gsl_blas_dnrm2(fxplus);

      if (fxplussize < (1-lambda/2.0)*fxsize) {
        gsl_vector_memcpy(x,xplus);
        f(x,fx);
        break;
      } else {
        lambda /=2.0;
      }

    } while(lambda > 1.0/64);


  stepsdone++;
  //printf("Steps done %i\n", stepsdone);
} while( gsl_blas_dnrm2(fx)>epsilon && stepsdone < maxsteps);

  gsl_vector_memcpy(xstart,x);
  fprintf(stderr, "Newton steps done %i\n", stepsdone);

  gsl_matrix_free(J);
  gsl_matrix_free(R);
  gsl_vector_free(x);
  gsl_vector_free(fx);
  gsl_vector_free(mfx);
  gsl_vector_free(xplus);
  gsl_vector_free(fxplus);
  gsl_vector_free(deltax);
}

void newtonWJ(
              void f(gsl_vector * X, gsl_vector * FX, gsl_matrix * JACOBIAN),
              gsl_vector * xstart,
              double epsilon) {
  int n = xstart->size;
  gsl_vector * x = gsl_vector_calloc(n);
  gsl_vector * xplus = gsl_vector_calloc(n);
  gsl_vector * deltax = gsl_vector_calloc(n);
  gsl_vector * fx = gsl_vector_calloc(n);
  gsl_vector * mfx = gsl_vector_calloc(n);
  gsl_matrix * J = gsl_matrix_calloc(n,n);
  gsl_matrix * R = gsl_matrix_calloc(n,n);

  double fxsize;
  double fxplussize;
  double lambda;
  int stepsdone = 0;
  int maxsteps = 1000000;

  gsl_vector_memcpy(x, xstart);
  do {
    f(x,fx,J);

    gsl_vector_memcpy(mfx, fx);
    gsl_vector_scale(mfx,-1);

    qrdecomp(J,R);
    qrbacksub(J,R,mfx,deltax);

    lambda = 1.0;
    do {
      fxsize = gsl_blas_dnrm2(fx);

      gsl_vector_memcpy(xplus,x);
      gsl_vector_scale(deltax,lambda);
      gsl_vector_add(xplus,deltax);
      gsl_vector_scale(deltax, 1.0/lambda);
      f(xplus,fx,J);

      fxplussize = gsl_blas_dnrm2(fx);

      if (fxplussize < (1-lambda/2.0)*fxsize) {
        gsl_vector_memcpy(x,xplus);
        break;
      } else {
        lambda /=2.0;
      }

    } while(lambda > 1.0/64);


  stepsdone++;
  //printf("Steps done %i\n", stepsdone);
} while( gsl_blas_dnrm2(fx)>epsilon && stepsdone < maxsteps);
  gsl_vector_memcpy(xstart,x);
  fprintf(stderr, "Newton w. Jacobian steps done %i\n", stepsdone);

  gsl_vector_free(x);
  gsl_vector_free(fx);
  gsl_vector_free(xplus);
  gsl_vector_free(deltax);
  gsl_vector_free(mfx);
  gsl_matrix_free(J);
  gsl_matrix_free(R);
}
