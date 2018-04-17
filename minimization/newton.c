#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

void qrdecomp(gsl_matrix * A, gsl_matrix * R);
void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x);

int newton(double f(gsl_vector * X, gsl_vector * DF, gsl_matrix * HESSIAN),
            gsl_vector * xstart,
            double eps) {
  int n = xstart->size;
  int steps = 0;
  int maxsteps = 1e3;
  double alpha = 1e-6;
  gsl_vector * dx = gsl_vector_calloc(n);
  gsl_vector * x = gsl_vector_calloc(n);
  gsl_vector * xnew = gsl_vector_calloc(n);
  gsl_vector * df = gsl_vector_calloc(n);
  gsl_vector * dfnew = gsl_vector_calloc(n);
  gsl_matrix * H = gsl_matrix_calloc(n,n);
  gsl_matrix * R = gsl_matrix_calloc(n,n);
  double lambda;
  double dfsize;
  double fval;
  double newfval;
  double dotprod;

  gsl_vector_memcpy(x,xstart);
  fval=f(x, df, H);

  do {

    gsl_vector_scale(df,-1);
    qrdecomp(H,R);
    qrbacksub(H,R,df,dx);
    gsl_vector_scale(df,-1);

    lambda = 1.0;
    do {
      gsl_vector_memcpy(xnew,x);
      gsl_vector_scale(dx,lambda);
      gsl_vector_add(xnew,dx);

      newfval = f(xnew, dfnew, R);

      gsl_blas_ddot(dx,df, &dotprod);

      if (newfval<fval+alpha*dotprod) {
        gsl_vector_memcpy(x,xnew);
        f(x,df,H);
        break;
      } else {
        gsl_vector_scale(dx,1/lambda);
        lambda /=2;
      }
    } while(lambda > 1.0/512);
    dfsize = gsl_blas_dnrm2(df);
    steps++;
    //fprintf(stderr, "dfsize = %g, eps = %g \n", dfsize, eps);
  } while(dfsize > eps && steps < maxsteps);

  gsl_vector_memcpy(xstart,x);

  gsl_vector_free(dx);
  gsl_vector_free(x);
  gsl_vector_free(xnew);
  gsl_vector_free(df);
  gsl_vector_free(dfnew);
  gsl_matrix_free(H);
  gsl_matrix_free(R);
  return steps;
}
