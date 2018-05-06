#include<math.h>
#include<assert.h>
#include<stdio.h>

int maxrecs = 1000000;


double integrator24(double f(double), double a, double b, double acc, double eps, double f2, double f3, double *interr,int nrecs){
  assert(nrecs<maxrecs);
  double f1 = f(a+(b-a)/6);
  double f4 = f(a+5*(b-a)/6);
  double goodQ = (2*f1+f2+f3+2*f4)/6 * (b-a);
  double badQ = (f1+f2+f3+f4)/4 * (b-a);



  double tol = acc+eps*fabs(goodQ);
  double error = fabs(goodQ - badQ);
  if (error < tol) {
    *interr = error;
    return goodQ;
  } else {
    double interr1;
    double interr2;

    double Q1 = integrator24(f, a, (a+b)/2, acc/sqrt(2.0), eps, f1, f2, &interr1,nrecs+1);
    double Q2 = integrator24(f, (a+b)/2, b, acc/sqrt(2.0), eps, f3, f4, &interr2,nrecs+1);

    *interr = sqrt(interr1*interr1 + interr2*interr2);
    return Q1+Q2;
  }
}

double integrator(double f(double), double a, double b, double acc, double eps, double *err){
  int nrecs=0;

  double f2 = f(a+2*(b-a)/6);
  double f3 = f(a+4*(b-a)/6);

  double Q =  integrator24(f, a, b, acc, eps, f2, f3, err, nrecs);
  printf("Q=%g\n", Q);
	//*err = my_estimate_of_the_error; /* should be smaller than acc+|Q|*eps */
	return Q;
}

double infintegrator(double f(double), double a, double b, double acc, double eps, double *err){
  int ainf = isinf(a);
  int binf = isinf(b);
  double Q;

  if (ainf ==0 && binf ==0) {
    //No infinities, business as usual
    Q = integrator( f, a, b, acc, eps, err);
  } else if (ainf != 0 && binf != 0) {
    //Infinity at both ends, gotta fix
    double Dinf(double t){return f(t/(1-t*t))*(1+t*t)/((1-t*t)*(1-t*t));}
    Q = integrator( Dinf, -1, 1, acc, eps, err);
  } else if (ainf == 0 && binf != 0){
    //Infinity at top of interval, fix
    double Uinf(double t){return f(a+(1-t)/t)/(t*t);}
    Q = integrator( Uinf, 0, 1, acc, eps, err);
  } else if (ainf != 0 && binf ==0) {
    //Infinity at bottom of interval, fix
    double Linf(double t){return f(b-(1-t)/t)/(t*t);}
    Q = integrator( Linf, 0, 1, acc, eps, err);
  } else {
    printf("Error in infinity function!\n");
    Q = INFINITY;
  }
	return Q;
}

double ClenshawCurtis(double f(double), double a, double b, double acc, double eps, double *err){
  double g(double t){return f( (a+b)/2+(a-b)/2*cos(t) )*sin(t)*(b-a)/2;}
  return integrator(g,0,M_PI,acc,eps,err);
}
