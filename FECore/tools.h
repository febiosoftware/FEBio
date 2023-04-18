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
#include "fecore_api.h"
#include <vector>

FECORE_API void linmin(double* p, double* xi, int n, double* fret, double(*fnc)(double[]));
FECORE_API void powell(double* p, double* xi, int n, double ftol, int* iter, double* fret, double(*fnc)(double[]));
FECORE_API double brent(double ax, double bx, double cx, double(*f)(double), double tol, double* xmin);
FECORE_API void mnbrak(double* ax, double* bx, double* cx, double* fa, double* fb, double* fc, double(*fnc)(double));
FECORE_API double golden(double ax, double bx, double cx, double(*f)(double), double tol, double* xmin);
FECORE_API double zbrent(double f(double, void*), double x1, double x2, double tol, void* data);
FECORE_API bool zbrac(double f(double, void*), double& x1, double& x2, void* data);
FECORE_API void solve_3x3(double A[3][3], double b[3], double x[3]);

FECORE_API bool LinearRegression(const std::vector<std::pair<double, double> >& data, std::pair<double, double>& res);
FECORE_API bool NonlinearRegression(const std::vector<std::pair<double, double> >& data, std::vector<double>& res, int func);
