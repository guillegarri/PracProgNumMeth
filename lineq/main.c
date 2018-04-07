#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdlib.h>

void qrdecomp(gsl_matrix * A, gsl_matrix * R);
void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x);
void inverse(gsl_matrix * A, gsl_matrix * B);


double myrandom(){
  return ((double)rand()/RAND_MAX);
}

void printM(gsl_matrix * A) {
  for (int i = 0; i < A->size1; i++) {
    for (int j = 0; j < A->size2; j++) {
      printf("%g ",gsl_matrix_get(A,i,j) );
    }
    printf("\n");
  }
}

int main(int argc, char const *argv[]) {
  // QR decomposition:
  int n=5;
  int m=4;

time_t t;
srand((unsigned) time(&t));

  gsl_matrix * A = gsl_matrix_calloc(n,m);
  gsl_matrix * R = gsl_matrix_calloc(m,m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      gsl_matrix_set(A,i,j,myrandom());
    }
  }

  printf("Random matrix A:\n");
  printM(A);

  printf("Performing QR decomposition\n");
  qrdecomp(A,R);
  printf("Resulting matrix Q:\n");
  printM(A);
  printf("Resulting matrix R:\n");
  printM(R);

  gsl_matrix * QTQ = gsl_matrix_calloc(m,m);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, A,A,0.0, QTQ );
  printf("Q^TQ should be 1, it is:\n");
  printM(QTQ);

  gsl_matrix * QR = gsl_matrix_calloc(n,m);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, R, 0.0, QR);
  printf("QR should be equal to A, it is:\n");
  printM(QR);

  // Linear system:

  m=5;

  gsl_matrix * SYS = gsl_matrix_calloc(m,m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      gsl_matrix_set(SYS,i,j,myrandom());
    }
  }

  printf("Square random matrix:\n");
  printM(SYS);

  gsl_matrix * QQ = gsl_matrix_alloc(m,m);
  gsl_matrix_memcpy(QQ,SYS);

  gsl_vector * b = gsl_vector_alloc(m);
  for (int i = 0; i < m; i++) {
    gsl_vector_set(b,i,myrandom());
  }
  printf("Random right hand side b:\n");
  gsl_vector_fprintf(stdout, b, "%7.3f");

  gsl_matrix * RR = gsl_matrix_calloc(m,m);

  qrdecomp(QQ,RR);

  gsl_vector * x = gsl_vector_alloc(m);
  qrbacksub(QQ,RR,b,x);

  printf("The solution x is:\n");
  gsl_vector_fprintf(stdout, x, "%7.3f");

  printf("We should have Ax=b, we have Ax=\n");
  gsl_blas_dgemv(CblasNoTrans, 1.0, SYS, x, 0.0, b);
  gsl_vector_fprintf(stdout, b, "%7.3f");


  // Matrix inverse
  
  gsl_matrix * Ainv = gsl_matrix_calloc(m,m);
  gsl_matrix_memcpy(QQ,SYS);
  inverse(QQ,Ainv);

  printf("The inverse matrix is:\n");
  printM(Ainv);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, SYS, Ainv, 0.0, RR);
  printf("AA^-1 should be 1, it is:\n");
  printM(RR);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ainv, SYS, 0.0, RR);
  printf("A^-1 A should be 1, it is:\n");
  printM(RR);


  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_matrix_free(QTQ);
  gsl_matrix_free(QR);
  gsl_matrix_free(SYS);
  gsl_matrix_free(QQ);
  gsl_vector_free(b);
  gsl_matrix_free(RR);
  gsl_vector_free(x);
  gsl_matrix_free(Ainv);
  return 0;
}
