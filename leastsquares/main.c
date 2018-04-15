#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include<math.h>
#include <stdio.h>

double funs(int i, double x){
   switch(i){
   case 0: return log(x); break;
   case 1: return 1.0;   break;
   case 2: return x;     break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}

void printM(gsl_matrix * A) {
  for (int i = 0; i < A->size1; i++) {
    for (int j = 0; j < A->size2; j++) {
      printf("%g ",gsl_matrix_get(A,i,j) );
    }
    printf("\n");
  }
}

void qrdecomp(gsl_matrix * A, gsl_matrix * R);
void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x);
void inverse(gsl_matrix * A, gsl_matrix * B);

int main(int argc, char const *argv[]) {
  int n, m;
  n = 9;
  m = 3;
  gsl_matrix * A = gsl_matrix_calloc(n,m);
  gsl_matrix * Q = gsl_matrix_calloc(n,m);
  gsl_vector * c = gsl_vector_calloc(m);
  gsl_vector * b = gsl_vector_calloc(n);
  double x[9] = {0.1, 1.33, 2.55, 3.78, 5.0, 6.22, 7.45, 8.68, 9.9};
  double y[9] = {-15.3, 0.32, 2.45, 2.75, 2.27, 1.35, 0.157, -1.23, -2.75};
  double dy[9] = {1.04, 0.594, 0.983, 0.998, 1.11, 0.398, 0.535, 0.968, 0.478};

  for (int i = 0; i < n; i++) {
    fprintf(stderr, "%g %g %g\n",x[i],y[i],dy[i]);
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      gsl_matrix_set(A,i,j,funs(j,x[i])/dy[i]);
    }
  }

  gsl_matrix_memcpy(Q, A);

  for (int i = 0; i < n; i++) {
    gsl_vector_set(b,i,y[i]/dy[i]);
  }

  gsl_matrix * R = gsl_matrix_calloc(m,m);

  qrdecomp(Q,R);
  qrbacksub(Q,R,b,c);

  gsl_matrix * Rinv = gsl_matrix_calloc(m,m);
  inverse(R,Rinv);

  gsl_matrix * S = gsl_matrix_calloc(m,m);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Rinv, Rinv, 0.0, S);

  printf("Resulting fit is (%g +- %g)*log(x)+(%g+-%g)+(%g+-%g)*x\n",
          gsl_vector_get(c,0), pow(gsl_matrix_get(S,0,0),0.5),
          gsl_vector_get(c,1), pow(gsl_matrix_get(S,1,1),0.5),
          gsl_vector_get(c,2), pow(gsl_matrix_get(S,2,2),0.5));

  double X, Y, Ymin, Ymax;
  fprintf(stderr, "\n\n");
  for (int i = 0; i < 1000; i++) {
    X = i/100.0;
    Y = gsl_vector_get(c,0)*log(X) + gsl_vector_get(c,1) + gsl_vector_get(c,2)*X;
    Ymin =  (gsl_vector_get(c,0)-pow(gsl_matrix_get(S,0,0),0.5))*log(X) +
            (gsl_vector_get(c,1)-pow(gsl_matrix_get(S,1,1),0.5)) +
            (gsl_vector_get(c,2)-pow(gsl_matrix_get(S,2,2),0.5))*X;
    Ymax =  (gsl_vector_get(c,0)+pow(gsl_matrix_get(S,0,0),0.5))*log(X) +
            (gsl_vector_get(c,1)+pow(gsl_matrix_get(S,1,1),0.5)) +
            (gsl_vector_get(c,2)+pow(gsl_matrix_get(S,2,2),0.5))*X;
    fprintf(stderr, "%g %g %g %g \n",X,Y,Ymin,Ymax);

  }






  gsl_matrix_free(A);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);
  gsl_matrix_free(Rinv);
  gsl_matrix_free(S);
  gsl_vector_free(c);
  gsl_vector_free(b);
  return 0;
}
