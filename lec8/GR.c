#include<stdio.h>
#include<math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>



int orbit_ode(double x, const double y[], double dydx[], void* params){
  double epsilon = *(double *) params;
  dydx[0] = y[1];
  dydx[1] = 1 - y[0] + epsilon*y[0]*y[0];
  return GSL_SUCCESS;
}



double my_orbit(double x, double epsilon, double val){
  gsl_odeiv2_system ucalc;
    ucalc.function=orbit_ode;
    ucalc.jacobian=NULL;
	  ucalc.dimension=2;
    ucalc.params=(void *) &epsilon;

  double hstart=copysign(0.1,x);
  double acc=1e-6;
  double eps=1e-6;
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&ucalc,gsl_odeiv2_step_rkf45,hstart,acc,eps);

  double t=0;
  double y[] = {1, val};
  gsl_odeiv2_driver_apply(driver,&t,x,y);

	gsl_odeiv2_driver_free(driver);
	return y[0];
}

int main() {
  for (int x = 0; x < 5*629; x++) {
    printf("%g %g \n", x/100.0, my_orbit(x/100.0, 0, 0) );
  }
  printf("\n");
  printf("\n");
  for (int x = 0; x < 5*629; x++) {
    printf("%g %g \n", x/100.0, my_orbit(x/100.0, 0, -0.5) );
  }
  printf("\n");
  printf("\n");
  for (int x = 0; x < 5*629; x++) {
    printf("%g %g \n", x/100.0, my_orbit(x/100.0, 0.01, -0.5) );
  }

  return 0;
}
