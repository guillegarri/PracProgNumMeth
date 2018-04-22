#include <math.h>
#include <stdio.h>

int main() {
  for (int i = 0; i < 200*M_PI; i++) {
    printf("%g %g %g\n",i/100.0, sin(i/100.0),cos(i/100.0));
  }
  return 0;
}
