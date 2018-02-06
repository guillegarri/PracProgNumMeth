#include "stdio.h"
#include "math.h"
#include "complex.h"
#include "limits.h"
#include "float.h"

int main() {
  double x = 1;
  while(1+x!=1){x/=2;}
  x*=2;
  printf("For while double epsilon\t = %g\n", x);

  float y = 1;
  while(1+y!=1){y/=2;}
  y*=2;
  printf("For while float epsilon \t = %f\n", y);

  long double z = 1;
  while(1+z!=1){z/=2;}
  z*=2;
  printf("For while long double epsilon\t = %Lg\n", z);

  double a =1;
  for (a; 1+a!=1; a/=2) {}
  a*=2;
  printf("For for double epsilon\t\t = %g\n", a);

  float b =1;
  for (b; 1+b!=1; b/=2) {}
  b*=2;
  printf("For for float epsilon\t\t = %f\n", b);

  long double c =1;
  for (c; 1+c!=1; c/=2) {}
  c*=2;
  printf("For for long double epsilon\t = %Lg\n", c);

  double X = 1;
  do {X/=2;
  } while(1+X!=1);
  X*=2;
  printf("For do while double epsilon\t = %g\n", X);

  float Y = 1;
  do {Y/=2;
  } while(1+Y!=1);
  Y*=2;
  printf("For do while float epsilon\t = %f\n", Y);

  long double Z = 1;
  do {Z/=2;
  } while(1+Z!=1);
  Z*=2;
  printf("For do while long double epsilon = %Lg\n", Z);

  printf("Machine double epsilon \t\t = %g\n", DBL_EPSILON);
  printf("Machine float epsilon \t\t = %f\n", FLT_EPSILON);
  printf("Machine long double epsilon\t = %Lg\n", LDBL_EPSILON);
  return 0;
}
