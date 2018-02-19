#include<stdio.h>
#include"nvector.h"
#include"stdlib.h"

nvector* nvector_alloc(int n){
  nvector* v = malloc(sizeof(nvector));
  (*v).size = n;
  (*v).data = malloc(n*sizeof(double));
  if( v==NULL ) fprintf(stderr,"error in nvector_alloc\n");
  return v;
}

void nvector_free(nvector* v){ free(v->data); free(v);}

void nvector_set(nvector* v, int i, double value){ (*v).data[i]=value; }

double nvector_get(nvector* v, int i){return (*v).data[i]; }

double nvector_dot_product (nvector* u, nvector* v){
  double sum = 0;
  for (int i = 0; i < (*v).size ; i++) {
    sum += (*v).data[i] * (*u).data[i];
  }
  return sum;
}

void nvector_print(char* s, nvector* v) {
  printf("%s [", s);
  for (int i = 0; i < (*v).size -1 ; i++) {
    printf("%g, ",(*v).data[i]);
  }
  printf("%g] \n", (*v).data[(*v).size -1]);
}

void nvector_set_zero(nvector* v) {
  for (int i = 0; i < (*v).size ; i++) {
    (*v).data[i] = 0;
  }
}

int nvector_equal(nvector* a, nvector* b){
  if ((*a).size == (*b).size) {
    for (int i = 0; i < (*a).size ; i++) {
      if ((*a).data[i] == (*b).data[i]) {
        continue;
      } else {
        return 0;
      }
    }
    return 1;
  } else {
    return 0;
  }
}

void nvector_add(nvector* a, nvector* b){
  if ((*a).size == (*b).size) {
    for (int i = 0; i < (*a).size ; i++) {
      (*a).data[i] += (*b).data[i];
    }
  } else {
    printf("Error, vectors do not have same size.\n");
  }
}

void nvector_sub(nvector* a, nvector* b){
  if ((*a).size == (*b).size) {
    for (int i = 0; i < (*a).size ; i++) {
      (*a).data[i] -= (*b).data[i];
    }
  } else {
    printf("Error, vectors do not have same size.\n");
  }
}

void nvector_scale(nvector* a, double x){
    for (int i = 0; i < (*a).size ; i++) {
      (*a).data[i] *= x;
    }
}
