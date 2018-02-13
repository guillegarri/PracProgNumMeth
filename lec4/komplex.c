#include<stdio.h>
#include"komplex.h"
#include "math.h"

void komplex_print (char *s, komplex a) {
	printf ("%s (%g,%g)\n", s, a.re, a.im);
}

komplex komplex_new (double x, double y) {
	komplex z = { x, y };
	return z;
}

void komplex_set (komplex* z, double x, double y) {
	(*z).re = x;
	(*z).im = y;
}

komplex komplex_add (komplex a, komplex b) {
	komplex result = { a.re + b.re , a.im + b.im };
	return result;
}

komplex komplex_sub (komplex a, komplex b){
  komplex result = {a.re - b.re , a.im - b.im};
  return result;
}

/*
int komplex_equal (komplex a, komplex b) {
  if (a.re != b.re) {
    return 0;
  } else if (a.im != b.im) {
    return 0;
  } else {
    return 1;
  }
}*/

int komplex_equal(komplex a, komplex b, double tau, double epsilon){
  if (fabs(a.re-b.re)<tau && fabs(a.im-b.im)<tau) {
    return 1;
  } else if (fabs(a.re-b.re)/(fabs(a.re)+fabs(b.re)) < epsilon/2.0 && fabs(a.im-b.im)/(fabs(a.im)+fabs(b.im)) < epsilon/2.0) {
    return 1;
  } else {
    return 0;
  }

}

komplex komplex_mul (komplex a, komplex b) {
  komplex result = {a.re * b.re - a.im * b.im , a.im * b.re + a.re * b.im };
  return result;
}

komplex komplex_div (komplex a, komplex b) {
  komplex result = {(a.re * b.re + a.im * b.im)/(b.re * b.re + b.im * b.im), (a.im * b.re - a.re * b.im)/(b.re * b.re + b.im * b.im)};
  return result;
}

komplex komplex_conjugate (komplex z) {
  komplex result = {z.re, -z.im};
  return result;
}

double komplex_abs (komplex z) {
  double result = sqrt(z.re * z.re + z.im * z.im);
  return result;
}

komplex komplex_exp (komplex z) {
  komplex result = {cos(z.im) * exp(z.re) , sin(z.im) * exp(z.re)};
  return result;
}

komplex komplex_timesi (komplex z) {
  komplex result = {-z.im , z.re};
  return result;
}

komplex komplex_sin (komplex z) {
  komplex result = {sin(z.re) * cosh(z.im), cos(z.re) * sinh(z.im)};
  return result;
}

komplex komplex_cos (komplex z) {
  komplex result = {cos(z.re) * cosh(z.im) , - sin(z.re) * sinh(z.im)};
  return result;
}

komplex komplex_sqrt ( komplex z) {
  double phi = atan(z.im / z.re);
  komplex result = {sqrt(komplex_abs(z)) * cos(phi/2) , sqrt(komplex_abs(z)) * sin(phi/2) };
  return result;
}
