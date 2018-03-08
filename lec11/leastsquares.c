#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multimin.h>

typedef struct {int n; double *t,*y,*e;} experimental_data;

double function_to_minimize (const gsl_vector *x, void *params) {
	double  A = gsl_vector_get(x,0);
	double  T = gsl_vector_get(x,1);
	double  B = gsl_vector_get(x,2);
	experimental_data *p = (experimental_data*) params;
	int     n = p->n;
	double *t = p->t;
	double *y = p->y;
	double *e = p->e;
	double sum=0;
  double f(double t) {return A*exp(-(t)/T) + B;}
//	#define f(t) A*exp(-(t)/T) + B
	for(int i=0;i<n;i++) sum += pow( (f(t[i]) - y[i] )/e[i] ,2);
	return sum;
}

int main() {
  double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
  double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
  double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
  int n = sizeof(t)/sizeof(t[0]);

  experimental_data data;
  data.n=n;
  data.t=t;
  data.y=y;
  data.e=e;

  const int dim=3;
  gsl_multimin_fminimizer *M
	   =gsl_multimin_fminimizer_alloc
		   (gsl_multimin_fminimizer_nmsimplex2,dim);

   gsl_multimin_function F;
   F.f=function_to_minimize;
   F.n=dim;
   F.params=(void*)&data;

   gsl_vector* start = gsl_vector_alloc(dim);
   gsl_vector_set(start,0,4);
   gsl_vector_set(start,1,3);
   gsl_vector_set(start,2,1);

   gsl_vector* step = gsl_vector_alloc(dim);
   gsl_vector_set(step,0,0.1);
   gsl_vector_set(step,1,0.1);
   gsl_vector_set(step,2,0.1);

   gsl_multimin_fminimizer_set(M,&F,start,step);

   int iter =0;

   do{
   gsl_multimin_fminimizer_iterate(M);
   iter++;
//   double A=gsl_vector_get(M->x,0);
//   double T=gsl_vector_get(M->x,1);
//   double B=gsl_vector_get(M->x,2);
   int status = gsl_multimin_test_size(M->size,1e-2);
   if(status==GSL_SUCCESS)break;
   }while(1);

   double A=gsl_vector_get(M->x,0);
   double T=gsl_vector_get(M->x,1);
   double B=gsl_vector_get(M->x,2);

   double f(double t) {return A*exp(-(t)/T) + B;}

   for(double x=t[0];x<t[n-1]+0.1;x+=0.05)printf("%g %g\n",x,f(x));
   fprintf(stderr,"A=%g T=%g B=%g iter=%i\n",A,T,B,iter);

   gsl_vector_free(start);
   gsl_vector_free(step);
   gsl_multimin_fminimizer_free(M);

  return 0;
}
