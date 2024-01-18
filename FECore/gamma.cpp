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
// adapted from Numerical Recipes in Fortran

bool gser(double& gamser, double a, double x, double& gln) {
    const int ITMAX = 100;
    const double EPS = 3e-7;
    gln = gammaln(a);
    if (x > 0) {
        double ap = a;
        double sum = 1./a;
        double del = sum;
        bool cnvgd = true;
        for (int n=1; n<=ITMAX; ++n) {
            ap += 1;
            del*= x/ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS) break;
            if (n == ITMAX) cnvgd = false;
        }
        if (!cnvgd) return false;
        gamser = sum*exp(-x+a*log(x)-gln);
    }
    else {
        gamser = 0;
        return false;
    }
	return true;
}

bool gcf(double& gammcf, double a, double x, double& gln) {
    const int ITMAX = 100;
    const double EPS = 3e-7, FPMIN = 10*std::numeric_limits<double>::epsilon();
    gln = gammaln(a);
    double b = x+1-a;
    double c = 1./FPMIN;
    double d = 1./b;
    double h = d;
    bool cnvgd = true;
    for (int i=1; i<=ITMAX; ++i) {
        double an = -i*(i-a);
        b += 2;
        d = an*d+b;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = b + an/c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1./d;
        double del = d*c;
        h *= del;
        if (fabs(del-1) < EPS) break;
        if (i == ITMAX) cnvgd = false;
    }
    if (!cnvgd) return false;
    gammcf = exp(-x+a*log(x)-gln)*h;
    return true;
}

double gamma_inc_P(double a, double x)
{
    // returns the lower incomplete Gamma function
    if ((x >= 0) && (a >= 0)) {
        if (x < a+1) {
            double gamser, gln;
            gser(gamser,a,x,gln);
            return gamser*exp(gln);
        }
        else {
            double gammcf, gln;
            gcf(gammcf, a, x, gln);
            return (1-gammcf)*exp(gln);
        }
    }
    return 0.0;
}

double gamma_inc_Q(double a, double x)
{
    // returns the upper incomplete Gamma function
    if ((x >= 0) && (a >= 0)) {
        if (x < a+1) {
            double gamser, gln;
            gser(gamser,a,x,gln);
            return (1-gamser)*exp(gln);
        }
        else {
            double gammcf, gln;
            gcf(gammcf, a, x, gln);
            return gammcf*exp(gln);
        }
    }
    return 0.0;
}
