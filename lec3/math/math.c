#include "stdio.h"
#include "math.h"
#include "complex.h"
int main(){
	double a = tgamma(5);
	printf("The gamma function of 5 is %g\n", a);
	double b = j1(0.5);
	printf("J1(0.5) is %g\n", b);
	complex c = csqrt(-2);
	printf("Square root of -2 is %g +i*%g\n", creal(c), cimag(c));
	complex d = cexp(I*1);
	printf("Complex exponential e to the power of i is %g +i*%g\n", creal(d), cimag(d));
  complex e = cexp(I*M_PI);
  printf("Complex exponential e to the power of i times pi is %g +i*%g\n", creal(e), cimag(e));
  complex f = cpow(I,M_E);
  printf("i to the power of e is %g +i*%g\n", creal(f), cimag(f));
  float x = 0.1111111111111111111111111111;
  double y = 0.1111111111111111111111111111;
  long double z = 0.1111111111111111111111111111L;
  printf("My float becomes \t %5.25f\nMy double becomes \t %5.25g\nMy long double becomes \t %5.25Lf\n",x,y,z );
}
