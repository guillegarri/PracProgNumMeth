#include<stdio.h>
#include<math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int logistic_ode(double x, const double y[], double dydx[], void* params){
  dydx[0] = y[0] * (1.0 - y[0]);
  return GSL_SUCCESS;
}

double my_logistic(double x){
  gsl_odeiv2_system sys;
	sys.function=logistic_ode;
	sys.jacobian=NULL;
	sys.dimension=1;
	sys.params=NULL;

  double hstart=copysign(0.1,x);
  double acc=1e-6;
  double eps=1e-6;
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,acc,eps);

  double t=0;
  double y[1] = {0.5};
  gsl_odeiv2_driver_apply(driver,&t,x,y);

	gsl_odeiv2_driver_free(driver);
	return y[0];
}

int main() {
  for (int x = 0; x < 300; x++) {
    printf("%g %g \n", x/100.0, my_logistic(x/100.0) );
  }
  printf("\n");
  printf("\n");
  for (int x = 0; x < 30; x++) {
    printf("%g %g \n", x/10.0, 1/(1+exp(-x/10.0)) );
  }
  return 0;
}
