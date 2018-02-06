#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

int equal(double a, double b, double tau, double epsilon);

int main(int argc, char *argv[]) {
  double a = atof(argv[1]);
  double b = atof(argv[2]);
  double tau = atof(argv[3]);
  double epsilon = atof(argv[4]);
  char yes[] = "Yes.";
  char no[] = "No.";
  printf("With tau =%g and epsilon = %g ",tau, epsilon);
  int isequal = equal(a, b, tau, epsilon);
  if (isequal==1) {
    printf("are %g and %g equal? %s\n",a, b, yes);
  } else if (isequal==0) {
    printf("are %g and %g equal? %s\n",a, b, no);
  } else {
    printf("Error\n");
  }
  return 0;
}
