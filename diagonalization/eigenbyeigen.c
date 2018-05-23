#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h>
#include <stdio.h>

int jacobi(gsl_matrix * A, gsl_vector * e, gsl_matrix * V);
int jacobiN_LtH(gsl_matrix * A, gsl_vector * e, gsl_matrix * V, int numvals);
int jacobiN_HtL(gsl_matrix * A, gsl_vector * e, gsl_matrix * V, int numvals);


void reset_rotations();
int get_rotations();

void printM(gsl_matrix * A) {
  for (int i = 0; i < A->size1; i++) {
    for (int j = 0; j < A->size2; j++) {
      printf("%+-5.3f ",gsl_matrix_get(A,i,j) );
    }
    printf("\n");
  }
}

double myrandom(){
  return ((double)rand()/RAND_MAX);
}

int main(int argc, char const *argv[]) {
  int n=7;
  int nvals = 1;

  time_t t;
  srand((unsigned) time(&t));

  gsl_matrix * A = gsl_matrix_calloc(n,n);
  gsl_matrix * Acopy = gsl_matrix_calloc(n,n);
  gsl_matrix * Acopy2 = gsl_matrix_calloc(n,n);
  gsl_matrix * Acopy3 = gsl_matrix_calloc(n,n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double Aij = myrandom();
      gsl_matrix_set(A, i, j, Aij);
      gsl_matrix_set(A, j, i, Aij);
    }
  }

  gsl_matrix_memcpy(Acopy, A);
  gsl_matrix_memcpy(Acopy2, A);
  gsl_matrix_memcpy(Acopy3, A);

  printf("My matrix A is\n");
  printM(A);

  gsl_matrix * V = gsl_matrix_calloc(n,n);
  gsl_matrix * Vcheck = gsl_matrix_calloc(n,n);

  gsl_vector * e = gsl_vector_alloc(n);
  gsl_vector * eH = gsl_vector_alloc(n);
  gsl_vector * echeck = gsl_vector_alloc(n);
  gsl_vector * eAll = gsl_vector_alloc(n);

  reset_rotations();

  int testsweeps = jacobi(A, echeck, Vcheck);
  int testrots = get_rotations();
  reset_rotations();

  int numsweeps = jacobiN_LtH(Acopy, e, V, nvals);
  int rots = get_rotations();

  reset_rotations();

  int Hsweeps = jacobiN_HtL(Acopy2, eH, V, nvals);
  int rotsH = get_rotations();

  reset_rotations();

  int sweepsall = jacobiN_HtL(Acopy3, eAll, V, n);
  int rotsall = get_rotations();

  printf("Lowest eigenvalue\n");
  gsl_vector_fprintf(stdout,e,"%7.3f");

  printf("Highest eigenvalue\n");
  gsl_vector_fprintf(stdout,eH,"%7.3f");

  printf("All eigenvalues value by value\n");
  gsl_vector_fprintf(stdout,eAll,"%7.3f");

  printf("All eigenvalues with cyclic sweeps\n");
  gsl_vector_fprintf(stdout,echeck,"%7.3f");

  printf("Rotations to get lowest eigenvalue: %i, sweeps: %i\n", rots, numsweeps);
  printf("Rotations to get highest eigenvalue: %i, sweeps: %i\n", rotsH, Hsweeps);
  printf("Rotations to get all eigenvalues (value by value): %i, sweeps: %i\n", rotsall, sweepsall);
  printf("Rotations to get all eigenvalues (cyclic): %i, sweeps: %i\n", testrots, testsweeps);


  gsl_vector_free(e);
  gsl_vector_free(eH);
  gsl_vector_free(eAll);
  gsl_matrix_free(A);
  gsl_matrix_free(Acopy);
  gsl_matrix_free(Acopy2);
  gsl_matrix_free(Acopy3);
  gsl_matrix_free(V);
  return 0;
}
