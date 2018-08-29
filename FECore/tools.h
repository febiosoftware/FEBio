#pragma once
#include "fecore_api.h"

FECORE_API void linmin(double* p, double* xi, int n, double* fret, double(*fnc)(double[]));
FECORE_API void powell(double* p, double* xi, int n, double ftol, int* iter, double* fret, double(*fnc)(double[]));
FECORE_API double brent(double ax, double bx, double cx, double(*f)(double), double tol, double* xmin);
FECORE_API void mnbrak(double* ax, double* bx, double* cx, double* fa, double* fb, double* fc, double(*fnc)(double));
FECORE_API double golden(double ax, double bx, double cx, double(*f)(double), double tol, double* xmin);
FECORE_API double zbrent(double f(double, void*), double x1, double x2, double tol, void* data);
FECORE_API bool zbrac(double f(double, void*), double& x1, double& x2, void* data);
FECORE_API void solve_3x3(double A[3][3], double b[3], double x[3]);
