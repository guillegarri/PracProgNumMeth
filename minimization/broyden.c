#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

int broyden(double f(gsl_vector * X, gsl_vector * DF),
            gsl_vector * xstart,
            double eps){
  int n = xstart->size;
  int steps = 0;
  int maxsteps = 1e6;
  double alpha = 1e-6;
  double minlambda = 1.0/512;
  double minsize = 1e-6;
  gsl_vector * dx = gsl_vector_calloc(n);
  //gsl_vector * s = gsl_vector_calloc(n);
  gsl_vector * x = gsl_vector_calloc(n);
  gsl_vector * xnew = gsl_vector_calloc(n);
  gsl_vector * df = gsl_vector_calloc(n);
  gsl_vector * dfnew = gsl_vector_calloc(n);
  gsl_matrix * B = gsl_matrix_calloc(n,n);
  gsl_vector * y = gsl_vector_calloc(n);
  gsl_vector * u = gsl_vector_calloc(n);
  gsl_vector * By = gsl_vector_calloc(n);
  double lambda;
  double dfsize;
  double fval;
  double newfval;
  double dotprod;
  double dxdoty;
  gsl_vector_memcpy(x,xstart);
  fval=f(x, df);
  gsl_matrix_set_identity(B);


  do {
    gsl_vector_scale(df,-1);

    gsl_blas_dgemv(CblasNoTrans, 1, B, df, 0, dx);

    gsl_vector_scale(df,-1);

    //gsl_vector_memcpy(s,dx);
    lambda = 1.0;

    do {
      gsl_vector_memcpy(xnew,x);
      gsl_vector_scale(dx,lambda);
      gsl_vector_add(xnew,dx);

      newfval = f(xnew, dfnew);

      gsl_blas_ddot(dx,df, &dotprod);

      if (newfval<fval+alpha*dotprod) {
        break;
      }
      if (lambda < minlambda) {
        gsl_matrix_set_identity(B);
        break;
      }
      lambda /= 2;
    } while(1);

    gsl_vector_memcpy(y,dfnew);

    gsl_vector_sub(y,df);

    gsl_vector_memcpy(u,dx);

    gsl_blas_dgemv(CblasNoTrans, 1, B, y, 0, By);

    gsl_vector_sub(u,By);

    gsl_blas_ddot(dx,y, &dxdoty);
    if (fabs(dxdoty)>minsize) {
      gsl_blas_dger(1/dxdoty, u, dx, B);
    }
    gsl_vector_memcpy(x,xnew);
    fval=f(x,df);
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
  gsl_vector_free(By);
  gsl_vector_free(y);
  gsl_vector_free(u);
  gsl_matrix_free(B);
  return steps;
}
