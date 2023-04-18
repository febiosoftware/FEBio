/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "math.h"
#include "fecore_api.h"

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
FECORE_API double gamma(double x);
