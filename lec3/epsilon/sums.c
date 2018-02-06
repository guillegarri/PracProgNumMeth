#include "stdio.h"
#include "math.h"
#include "complex.h"
#include "limits.h"
#include "float.h"

int max = INT_MAX;

int main() {
  float sum_up_float = 0.0f;
  for (int N=1; N < max; N++) {
    sum_up_float+=1.0f/N;
  }
  sum_up_float+=1.0f/max;
  printf("Sum up float = %f\n", sum_up_float);

  float sum_down_float = 0.0f;
  for (int N=1; N < max; N++) {
    sum_down_float+=1.0f/(max-N);
  }
  printf("Sum down float = %f\n", sum_down_float);

  double sum_up_double = 0.0;
  for (int N=1; N < max; N++) {
    sum_up_double+=1.0/N;
  }
  sum_up_double+=1.0/max;
  printf("Sum up double = %3.8g\n", sum_up_double);

  float sum_down_double = 0.0;
  for (int N=1; N < max; N++) {
    sum_down_double+=1.0/(max-N);
  }
  printf("Sum down double = %3.8g\n", sum_down_double);
/*
  float sum_up_float = 0.0f;
  int N = 1;
  do {
    sum_up_float += 1.0f/N;
    N++;
  } while(N<max);
  printf("Sum up float = %f\n", sum_up_float);

  N=1;
  int diff = 0;
  float sum_down_float = 0.0f;
  do {
    diff = max - N;
    sum_down_float += 1.0f/diff;
    N++;
  } while(N<max);
  printf("Sum down float = %f\n", sum_down_float);
  if (equal(5.0,5.2,0.3,0.2)==1) {
    char yn[4] = {'Y','e','s','\0'};
  } else {
    char yn[3] = {'N','o','\0'};
  }*/


  return 0;
}
