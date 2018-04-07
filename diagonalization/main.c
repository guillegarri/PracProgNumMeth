#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h>
#include <stdio.h>

int jacobi(gsl_matrix * A, gsl_vector * e, gsl_matrix * V);

void printM(gsl_matrix * A) {
  for (int i = 0; i < A->size1; i++) {
    for (int j = 0; j < A->size2; j++) {
      printf("%+-5.3f ",gsl_matrix_get(A,i,j) );
    }
    printf("\n");
  }
}

void sandwich(gsl_matrix * A, gsl_matrix * B , gsl_matrix * BTAB) {
  int n = A->size1;
  gsl_matrix * AB = gsl_matrix_calloc(n,n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, AB); //AB
  gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, B, AB, 0.0, BTAB);
  gsl_matrix_free(AB);
}

double myrandom(){
  return ((double)rand()/RAND_MAX);
}

int main(int argc, char const *argv[]) {

int n;

if (argc > 1) {
  n = atoi(argv[1]);
} else {
  n = 7;
}


  time_t t;
  srand((unsigned) time(&t));

  gsl_matrix * A = gsl_matrix_calloc(n,n);
  gsl_matrix * Acopy = gsl_matrix_calloc(n,n);


  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double Aij = myrandom();
      gsl_matrix_set(A, i, j, Aij);
      gsl_matrix_set(A, j, i, Aij);
    }
  }

  gsl_matrix_memcpy(Acopy, A);

  printf("My matrix A is\n");
  printM(A);

  gsl_matrix * V = gsl_matrix_calloc(n,n);

  gsl_vector * e = gsl_vector_alloc(n);

  int numsweeps = jacobi(A, e, V);
  printf("Diagonalization done in %i sweeps\n",numsweeps);
  printf("This yields eigenvalues\n");
  gsl_vector_fprintf(stdout, e, "%7.3f");

  printf("And a matrix of eigenvectors\n");
  printM(V);

  gsl_matrix * D = gsl_matrix_calloc(n,n);
  for (int i = 0; i < n; i++) {
    gsl_matrix_set(D,i,i,gsl_vector_get(e,i));
  }

  gsl_matrix * Dcomp = gsl_matrix_calloc(n,n);
  sandwich(Acopy, V, Dcomp);

  printf("V^TAV should yield\n");
  printM(D);
  printf("It yields\n");
  printM(Dcomp);




  return 0;
}
