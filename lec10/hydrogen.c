#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include "stdio.h"
#include<assert.h>

int hydrogen_ode(double r, const double y[], double dydx[], void* params){
  double epsilon = *(double *) params;
  dydx[0] = y[1];
  dydx[1] = -2*(epsilon*y[0]+y[0]/r);
  return GSL_SUCCESS;
}

int
print_state (int iter, gsl_multiroot_fsolver * s)
{
  fprintf (stderr,"iter = %3i x = % .10g f(x) = % .10g \n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->f, 0));
					return 0;
}

double my_hydrogen(double r, double epsilon){
  assert(r>=0);
  const double rmin = 1e-3;
  if (r<rmin) return r-r*r;

  gsl_odeiv2_system system;
  system.function = hydrogen_ode;
  system.jacobian = NULL;
  system.dimension = 2;
  system.params = (void*) &epsilon;

  double hstart=copysign(1e-3,r);
  double acc=1e-6;
  double eps=1e-6;
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new
  (&system,gsl_odeiv2_step_rkf45,hstart,acc,eps);

  double t=rmin, y[] = {t-t*t, 1-2*t};
  int status = gsl_odeiv2_driver_apply (driver, &t, r, y);
  if (status != GSL_SUCCESS) fprintf (stderr,"Fe: odeiv2 error: %d\n", status);

  gsl_odeiv2_driver_free (driver);
	return y[0];
}

int M_func(const gsl_vector* current_point, void* params, gsl_vector* M){
  double rmax = *(double *) params;
  double epsilon = gsl_vector_get(current_point,0);
  assert(epsilon<0);
  gsl_vector_set(M, 0, my_hydrogen(rmax, epsilon));
  return GSL_SUCCESS;
}

int main(int argc, char const *argv[]) {
  const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *solver;

	int dim = 1;
  double rmax = 8.0;

	gsl_multiroot_function F;
	F.f=M_func;
	F.n=dim;
  F.params = (void*) &rmax;

  int status;
  int iter=0;

  double eps_init = -14;
  gsl_vector *X = gsl_vector_alloc(dim);
  gsl_vector_set(X, 0, eps_init);

  T = gsl_multiroot_fsolver_hybrids;
  solver = gsl_multiroot_fsolver_alloc(T, dim);
  gsl_multiroot_fsolver_set(solver, &F, X);

  print_state(iter,solver);


  do {
		iter++;
		status = gsl_multiroot_fsolver_iterate(solver);
		print_state(iter,solver);
		if (status)
			break;

		status = gsl_multiroot_test_residual(solver->f, 1e-3);

	} while(status==GSL_CONTINUE && iter < 10000);

  double epsilon=gsl_vector_get(solver->x,0);

  printf("epsilon = %g\n", epsilon);
	//printf("status = %s\n", gsl_strerror(status));

  printf("\n\n");

  for(double r=0; r<=rmax; r+=rmax/64) printf("%g %g %g\n",r,my_hydrogen(r,epsilon),r*exp(-r));


  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(X);

  return 0;
}
