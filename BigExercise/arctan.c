#include<stdio.h>
#include<math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int arctan_ode(double x, const double y[], double dydx[], void* params){
  dydx[0]=1.0/(x*x +1);
  return GSL_SUCCESS;
}

double my_arctan(double x){
  gsl_odeiv2_system sys;
  sys.function = arctan_ode;
  sys.jacobian = NULL;
  sys.dimension = 1;
  sys.params = NULL;

  double startstep = copysign(0.1, x);
  double acc = 1e-6;
  double eps = 1e-6;
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new
    (&sys,gsl_odeiv2_step_rkf45,startstep,acc,eps);

    double startx = 0;
    double y[1] = {0};

    gsl_odeiv2_driver_apply(driver,&startx,x,y);

    gsl_odeiv2_driver_free(driver);
    return y[0];
}

int main(void) {
  for (double x = 0.0; x < 4*M_PI; x+=0.01) {
    printf("%g %.10g %.10g %.10g\n",x,my_arctan(x),atan(x), my_arctan(x)-atan(x));
  }
  return 0;
}
