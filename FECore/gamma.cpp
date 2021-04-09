/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "stdafx.h"
#include "gamma.h"
#include <limits>
#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

//---------------------------------------------------------------------------------
// evaluate natural logarithm of gamma function
double gammaln(double x)
{
    static double c[6] = {
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5};
    double z = x;
    double gln = x + 5.5;
    gln -= (x+0.5)*log(gln);
    double s = 1.000000000190015;
    for (int i=0; i<6; ++i) s += c[i]/++z;
    gln = -gln + log(2.5066282746310005*s/x);
    return gln;
}

//---------------------------------------------------------------------------------
// evaluate the inverse of the gamma function
double gammainv(double x) {
    if (x >= 0) return exp(-gammaln(x));
    else return sin(M_PI*x)*exp(gammaln(1-x))/M_PI;
}

//---------------------------------------------------------------------------------
// evaluate the incomplete Gamma function
double gamma_inc_P(double a, double x)
{
    const int MAXITER = 100;
    const double eps = 10*std::numeric_limits<double>::epsilon();
    bool convgd = false;
    int i = 0;
    double d = gammainv(a+1);
    double P = d;
    while (!convgd) {
        ++i;
        d *= x/(a+i);
        P += d;
        if ((fabs(d) < eps*fabs(P)) || (i > MAXITER)) convgd = true;
    }
    return P*exp(-x)*pow(x, a);
}

double gamma_inc_Q(double a, double x)
{
    return 1 - gamma_inc_P(a,x);
}
