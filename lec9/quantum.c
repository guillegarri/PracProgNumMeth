#include<gsl/gsl_integration.h>
#include<math.h>
#include<stdio.h>

double psisq(double x,  void* params){
  double alpha=*(double*)params;
  return exp(-alpha*x*x);
}

double psisqwH(double x, void* params){
  double alpha=*(double*)params;
  return (-alpha*alpha*x*x/2 +alpha/2 + x*x/2)*exp(-alpha*x*x);
}

double normsq(double alpha){
  gsl_function f;
  f.function=&psisq;
  f.params=(void*)&alpha;

  size_t limit=10;
  double acc=1e-6,eps=1e-6,result,err;
	gsl_integration_workspace* workspace=
		gsl_integration_workspace_alloc(limit);
	gsl_integration_qagi(&f,acc,eps,limit,workspace,&result,&err);

	gsl_integration_workspace_free(workspace);

  fprintf(stderr,"called normsq\n");

  return result;
}

double expH(double alpha){
  gsl_function f;
  f.function=&psisqwH;
  f.params=(void*)&alpha;

  size_t limit=10;
  double acc=1e-6,eps=1e-6,result,err;
	gsl_integration_workspace* workspace=
		gsl_integration_workspace_alloc(limit);
	gsl_integration_qagi(&f,acc,eps,limit,workspace,&result,&err);

	gsl_integration_workspace_free(workspace);

  fprintf(stderr,"called expH\n");

  return result;
}

int main(){
  for (double a = 0.1; a < 5.0; a+=0.1) {
    double E = expH(a)/normsq(a);
    printf("%g %g\n",a, E);
  }
  return 0;
}
