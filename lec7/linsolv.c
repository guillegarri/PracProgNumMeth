#include "math.h"
#include "stdio.h"
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

int main(int argc, char const *argv[]) {
  gsl_matrix * A = gsl_matrix_alloc(3, 3);
  gsl_matrix * Acop = gsl_matrix_alloc(3, 3);

  gsl_matrix_set(A, 0, 0, 6.13);
  gsl_matrix_set(A, 0, 1, -2.90);
  gsl_matrix_set(A, 0, 2, 5.86);
  gsl_matrix_set(A, 1, 0, 5.37);
  gsl_matrix_set(A, 1, 1, -6.31);
  gsl_matrix_set(A, 1, 2, -3.89);
  gsl_matrix_set(A, 2, 0, -4.36);
  gsl_matrix_set(A, 2, 1, 1.00);
  gsl_matrix_set(A, 2, 2, 0.19);

  gsl_matrix_set(Acop, 0, 0, 6.13);
  gsl_matrix_set(Acop, 0, 1, -2.90);
  gsl_matrix_set(Acop, 0, 2, 5.86);
  gsl_matrix_set(Acop, 1, 0, 5.37);
  gsl_matrix_set(Acop, 1, 1, -6.31);
  gsl_matrix_set(Acop, 1, 2, -3.89);
  gsl_matrix_set(Acop, 2, 0, -4.36);
  gsl_matrix_set(Acop, 2, 1, 1.00);
  gsl_matrix_set(Acop, 2, 2, 0.19);

  gsl_vector * y = gsl_vector_alloc(3);
  gsl_vector * x = gsl_vector_alloc(3);


  gsl_vector_set(y,0,6.23);
  gsl_vector_set(y,1,5.37);
  gsl_vector_set(y,2,2.29);

  printf("y =\n");
  gsl_vector_fprintf(stdout, y, "%g");

  gsl_linalg_HH_solve(A, y, x);
  printf("The solution to A*x=y is\n");
  gsl_vector_fprintf(stdout, x, "%g");

  gsl_blas_dgemv(CblasNoTrans, 1.0, Acop, x, 0.0, y);
  printf("A*x =\n");
  gsl_vector_fprintf(stdout, y, "%g");

  gsl_matrix_free(A);
  gsl_matrix_free(Acop);
  gsl_vector_free(y);
  gsl_vector_free(x);
  return 0;
}
