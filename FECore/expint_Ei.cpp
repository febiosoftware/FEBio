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
#include "expint_Ei.h"
#include <limits>
#include <math.h>

// This is our homemade function for evaluating the exponential integral Ei(x), valid for
// negative and positive values of x.

double expint_Ei(double x)
{
    const int MAXITER = 100;
    if (x == 0) return -INFINITY;
    const double eps = 10*std::numeric_limits<double>::epsilon();
    const double gamma = 0.5772156649015328606065120900824024310421;
    double ei = 0;

    // use series expansion if x is sufficiently small
    if (fabs(x) < fabs(log(eps))) {
        // use series expansion if x is sufficiently small
        ei = gamma;
        ei += (x > 0) ? log(x) : log(-x);
        double d = x;
        int i=0;
        bool convgd = false;
        while (!convgd) {
            ++i;
            ei += d;
            if ((fabs(d) < eps*fabs(ei)) || (i > MAXITER)) convgd = true;
            else d *= (i*x)/pow(i+1,2);
        }
    }
    else {
        // use asymptotic expansion for sufficiently large x
        ei = 1;
        double d = 1./x;
        int i=0;
        bool convgd = false;
        while (!convgd) {
            ++i;
            ei += d;
            if ((fabs(d) < eps*fabs(ei)) || (i > MAXITER)) convgd = true;
            else d *= (i+1)/x;
        }
        ei *= exp(x)/x;
    }
    return ei;
}
