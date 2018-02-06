#include "stdio.h"
#include "math.h"
#include "complex.h"
#include "limits.h"
#include "float.h"

int main(){
  int i=1; while(i-1<i){i--;}
  printf("For while loop, my min int \t = %i\n",i);

  int j=1;
  for ( j; j > j-1; j--){}
  printf("For for loop, my min int \t = %i\n",j);

  int k = 1;
  int z = 1;
  do {z=k;
    k--;
  } while(z>k);
  printf("For do while loop, my min int \t = %i\n",z);
  printf("My computer min int is \t \t = %i\n",INT_MIN);
}
