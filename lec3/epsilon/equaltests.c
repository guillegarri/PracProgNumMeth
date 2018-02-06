#include "stdio.h"
#include "math.h"

int equal(double a, double b, double tau, double epsilon);

int main() {
  double mytau = 0.3;
  double myepsilon = 0.05;
  char yes[] = "Yes.";
  char no[] = "No.";
  printf("My tau =%g and my epsilon = %g\n",mytau, myepsilon );
  int isequal = equal(5.0, 5.2, mytau, myepsilon);
  if (isequal==1) {
    printf("Are 5.0 and 5.2 equal? %s\n", yes);
  } else if (isequal==0) {
    printf("Are 5.0 and 5.2 equal? %s\n", no);
  } else {
    printf("Error\n");
  }


  return 0;
}
