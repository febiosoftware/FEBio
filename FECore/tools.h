#pragma once

void linmin(double* p, double* xi, int n, double* fret, double(*fnc)(double[]));
void powell(double* p, double* xi, int n, double ftol, int* iter, double* fret, double(*fnc)(double[]));
double brent(double ax, double bx, double cx, double(*f)(double), double tol, double* xmin);
void mnbrak(double* ax, double* bx, double* cx, double* fa, double* fb, double* fc, double(*fnc)(double));
double golden(double ax, double bx, double cx, double(*f)(double), double tol, double* xmin);
double zbrent(double f(double), double x1, double x2, double tol);