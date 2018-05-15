#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_multimin.h>

typedef struct {
	int n;
	double (*f)(double,double);
	gsl_vector* data;
} ann2d;

ann2d* ann2d_alloc(int n,double(*f)(double, double)){
	ann2d* nw=malloc(sizeof(ann2d));
	nw->n=n;
	nw->f=f;
	nw->data=gsl_vector_alloc(5*n);
	return nw;
}

void ann2d_free(ann2d* nw){
	gsl_vector_free(nw->data);
	free(nw);
}

double ann2d_feed(ann2d*nw,double x, double y){
	double val=0;
	for(int i=0;i<nw->n;i++){
		double a=gsl_vector_get(nw->data,0*nw->n+i);
		double b=gsl_vector_get(nw->data,1*nw->n+i);
		double c=gsl_vector_get(nw->data,2*nw->n+i);
		double d=gsl_vector_get(nw->data,3*nw->n+i);
		double w=gsl_vector_get(nw->data,4*nw->n+i);
		val+=nw->f((x-a)/b, (y-c)/d)*w;
	}
	return val;
}

void ann2d_train(ann2d * nw, gsl_vector * X, gsl_vector * Y, gsl_matrix * Z) {
	int n = nw->data->size;
	gsl_vector * v = gsl_vector_alloc(n);
	gsl_vector * startstep = gsl_vector_alloc(n);

#define DIFF(x) pow(x,2)
  double score(const gsl_vector * v, void * params){
    gsl_vector_memcpy(nw->data,v);

    double s=0;
    for (int i = 0; i < X->size; i++) {
			for (int j = 0; j < Y->size; j++) {
				double x = gsl_vector_get(X,i);
	      double y = gsl_vector_get(Y,j);
				double z = gsl_matrix_get(Z,i,j);
	      double fit = ann2d_feed(nw,x,y);
				//fprintf(stderr, "x = %g, y = %g, fit = %g\n",x,y,fit );
	      s += DIFF(z-fit);
			}
    }
    s /= X->size*Y->size;

		//fprintf(stderr, "score: %g\n", s);
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
			//fprintf(stderr, "%i iterations done, cannot be improved.\n", iters);
			break;
		}

		if (minimizer->size < 1e-5) {
			//fprintf(stderr, "converged\n");
			break;
		}
	} while(iters < 100000);

	gsl_vector_memcpy(nw->data, minimizer->x);

	gsl_vector_free(v);
	gsl_vector_free(startstep);
	gsl_multimin_fminimizer_free(minimizer);

}
