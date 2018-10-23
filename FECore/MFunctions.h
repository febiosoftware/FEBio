#pragma once
#include "math.h"

typedef double (*FUNCPTR) (double);
typedef double (*FUNC2PTR) (double, double);
typedef double (*FUNCNPTR) (double*, int);

// trigonometric functions
double csc(double x);
double sec(double x);
double cot(double x);
double sinc(double x);

// Bessel functions
#ifdef WIN32
double jn(double x, double y);
double yn(double x, double y);
#endif

// multi-variate functions
double fmax(double* x, int n);
double fmin(double* x, int n);
double avg(double* x, int n);

// special polynomials
double chebyshev(double f, double x);

// additional functios
double fl(double x);
double sgn(double x);
double fac(double x);	// factorials
double heaviside(double x);
double unit_step(double x);

// binomials
double binomial(double n, double r);

// gamma function
double gamma(double x);
