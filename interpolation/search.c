#include "stdio.h"

int search(int n, double * x, double z){
  //fprintf(stderr, "running search\n");
  int L = 0;
  int R = n-1;
  if (L>R) {
    printf("List must have length greater than 0, provided length is n=%i\n",n);
    return 1e10;
  }
  int iter=0;
  //fprintf(stderr, "doing search loop\n");
  do {


    int m = (L+R)/2;
    if (x[m] <= z && x[m+1] > z) {
    //  fprintf(stderr, "found m=%i\n",m);
      return m;
    } else if (x[m] < z) {
      L = m+1;
    } else if (x[m] > z) {
      R = m-1;
    }


    iter++;
    if (iter>n) {
      printf("search failed\n");
      return 1e10;
      break;
    }
    //fprintf(stderr, "did a loop, m=%i, x[m]=%g, z=%g\n",m,x[m],z);
  } while(1);
}
