#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void qrdecomp(gsl_matrix * A, gsl_matrix * R) {
  int m = A->size2;
  for (int i = 0; i < m; i++) {
    gsl_vector_view colAi = gsl_matrix_column(A,i);
    double collength = gsl_blas_dnrm2(&colAi.vector);
    gsl_matrix_set(R,i,i,collength);
    gsl_vector_scale(&colAi.vector, 1.0/collength);
    for (int j = i+1; j < m; j++) {
      gsl_vector_view colAj = gsl_matrix_column(A,j);
      double dotprod=0;
      gsl_blas_ddot(&colAi.vector, &colAj.vector, &dotprod);
      gsl_blas_daxpy(-dotprod, &colAi.vector, &colAj.vector);
      gsl_matrix_set(R,i,j,dotprod);
    }
  }
}

void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x) {
  int m=R->size1;
  gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, x);
  for (int i = m-1; i >=0; i--) {
    double s=0;
    for (int j = i+1; j < m; j++) {
      s += gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
    }
    gsl_vector_set(x,i,(gsl_vector_get(x,i)-s)/gsl_matrix_get(R,i,i));
  }
}

void inverse(gsl_matrix * A, gsl_matrix * B) {
  int n = A->size1;
  gsl_matrix * R = gsl_matrix_calloc(n,n);
  qrdecomp(A,R);
  gsl_vector * b = gsl_vector_calloc(n);
  gsl_vector * x = gsl_vector_calloc(n);
  for (int i = 0; i < n; i++) {
    gsl_vector_set(b,i,1.0);
    qrbacksub(A,R,b,x);
    gsl_vector_set(b,i,0.0);
    gsl_matrix_set_col(B,i,x);
  }

  gsl_matrix_free(R);
  gsl_vector_free(b);
  gsl_vector_free(x);
}
