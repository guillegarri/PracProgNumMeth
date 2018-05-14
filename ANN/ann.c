#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_multimin.h>

typedef struct {
	int n;
	double (*f)(double);
	gsl_vector* data;
	} ann;

ann* ann_alloc(int n,double(*f)(double)){
	ann* nw=malloc(sizeof(ann));
	nw->n=n;
	nw->f=f;
	nw->data=gsl_vector_alloc(3*n);
	return nw;
}

void ann_free(ann* nw){
	gsl_vector_free(nw->data);
	free(nw);
}

double ann_feed(ann*nw,double x){
	double val=0;
	for(int i=0;i<nw->n;i++){
		double a=gsl_vector_get(nw->data,0*nw->n+i);
		double b=gsl_vector_get(nw->data,1*nw->n+i);
		double w=gsl_vector_get(nw->data,2*nw->n+i);
		val+=nw->f((x-a)/b)*w;
	}
	return val;
}

void ann_train(ann * nw, gsl_vector * X, gsl_vector * Y) {
	int n = 3*nw->n;
	gsl_vector * v = gsl_vector_alloc(n);
	gsl_vector * startstep = gsl_vector_alloc(n);

#define DIFF(x) pow(x,2)
  double score(const gsl_vector * v, void * params){
    gsl_vector_memcpy(nw->data,v);

    double s=0;
    for (int i = 0; i < X->size; i++) {
      double x = gsl_vector_get(X,i);
      double y = gsl_vector_get(Y,i);
      double fit = ann_feed(nw,x);
			fprintf(stderr, "x = %g, y = %g, fit = %g\n",x,y,fit );
      s += DIFF(y-fit);
    }
//    s /= X->size;

		fprintf(stderr, "score: %g\n", s);
		gsl_vector_fprintf(stderr,v,"%8.3g");
    return s;
  }

	gsl_vector_memcpy(v, nw->data);

	gsl_multimin_function F = {.f=score, .n=n, .params=NULL};
	gsl_vector_set_all(startstep, 0.1);

  gsl_multimin_fminimizer * minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,n);
	gsl_multimin_fminimizer_set(minimizer, &F, v, startstep);

	int iters=0;
	//int status;
	int flag;
	do {
		iters++;
		flag = gsl_multimin_fminimizer_iterate(minimizer);
		if (flag!=0){
			fprintf(stderr, "%i iterations done, cannot be improved.\n", iters);
			break;
		}

		if (minimizer->size < 1e-4) {
			fprintf(stderr, "converged\n");
			break;
		}
	} while(iters < 100000);

	gsl_vector_memcpy(nw->data, minimizer->x);

	gsl_vector_free(v);
	gsl_vector_free(startstep);
	gsl_multimin_fminimizer_free(minimizer);

}
