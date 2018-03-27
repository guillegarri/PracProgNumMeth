#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"
typedef struct {int n; double * x, *y, *b, *c;} qspline;

qspline * qspline_alloc(int n, double * x, double * y);
double qspline_eval(qspline * s, double z);
double qspline_deriv(qspline * s, double z);
double qspline_int(qspline * s, double z);
void qspline_free(qspline * s);

int main() {
double x[] = {1, 2, 3, 4, 5};
double y1[] = {1, 1, 1, 1, 1};
double y2[] = {1, 2, 3, 4, 5};
double y3[] = {1, 2, 9, 16, 25};

qspline * s1 = qspline_alloc(5, x, y1);
qspline * s2 = qspline_alloc(5, x, y2);
qspline * s3 = qspline_alloc(5, x, y3);
printf("For the first interval \n");
for (int i = 0; i < 4; i++) {
  printf("Analytic b[%i]=0, program b[%i]=%g\n",i,i,s1->b[i]);
  printf("Analytic c[%i]=0, program c[%i]=%g\n",i,i,s1->c[i]);
}

printf("For the second interval \n");
for (int i = 0; i < 4; i++) {
  printf("Analytic b[%i]=1, program b[%i]=%g\n",i,i,s2->b[i]);
  printf("Analytic c[%i]=0, program c[%i]=%g\n",i,i,s2->c[i]);
}

double b3[4], c3[4];

c3[0]=0;
for (int i = 0; i < 3; i++) {
  c3[i+1] = ((y3[i+2]-y3[i+1])-(y3[i+1]-y3[i])-c3[i]);
}
c3[3] /= 2;

for (int i = 2; i >=0; i--) {
  c3[i] = ((y3[i+2]-y3[i+1])-(y3[i+1]-y3[i])-c3[i+1]);
}

for (int i = 0; i < 4; i++) {
  b3[i] = (y3[i+1]-y3[i])-c3[i];
}

printf("For the third interval \n");
for (int i = 0; i < 4; i++) {
  printf("Analytic b[%i]=0, Manual b[%i]=%g, program b[%i]=%g\n",i,i,b3[i],i,s3->b[i]);
  printf("Analytic c[%i]=1, Manual b[%i]=%g, program c[%i]=%g\n",i,i,c3[i],i,s3->c[i]);
}

qspline_free(s1);
qspline_free(s2);
qspline_free(s3);
  return 0;
}
