#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multimin.h>


double rosenbrock(const gsl_vector* X, void* params){
  double x = gsl_vector_get(X,0);
  double y = gsl_vector_get(X,1);
  return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}

int main(int argc, char const *argv[]) {
  const int dim=2;
  gsl_multimin_fminimizer *M
	=gsl_multimin_fminimizer_alloc
		(gsl_multimin_fminimizer_nmsimplex2,dim);

  gsl_multimin_function F;
  F.f=rosenbrock;
  F.n=dim;

  gsl_vector* start = gsl_vector_alloc(dim);
  gsl_vector_set(start,0,2);
  gsl_vector_set(start,1,2);

  gsl_vector* step = gsl_vector_alloc(dim);
  gsl_vector_set(step,0,0.1);
  gsl_vector_set(step,1,0.1);

  gsl_multimin_fminimizer_set(M,&F,start,step);


  int iter=0;
do{
  gsl_multimin_fminimizer_iterate(M);
  iter++;
  double x=gsl_vector_get(M->x,0);
  double y=gsl_vector_get(M->x,1);
  fprintf(stderr,"%g %g\n",x,y);
  int status = gsl_multimin_test_size(M->size,1e-2);
  if(status==GSL_SUCCESS)break;
}while(1);

double x = gsl_vector_get(M->x,0);
double y = gsl_vector_get(M->x,1);

gsl_vector_free(start);
gsl_vector_free(step);
gsl_multimin_fminimizer_free(M);

  printf("x_min=%g y_min=%g iterations=%i\n",x,y,iter );

  return 0;
}
