#include<gsl/gsl_integration.h>
#include<math.h>
#include<stdio.h>

double integrand(double x, void* params){
  return log(x)/sqrt(x);
}

int main(int argc, char const *argv[]) {
  gsl_function f;
  f.function=&integrand;
  f.params=NULL;

  size_t limit=10;
  double acc=1e-6,eps=1e-6,result,err;
	gsl_integration_workspace* workspace=
		gsl_integration_workspace_alloc(limit);
	gsl_integration_qags(&f,0,1,acc,eps,limit,workspace,&result,&err);

	gsl_integration_workspace_free(workspace);
  printf("The integral is %g with error %g\n", result, err);
  return 0;
}
