/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2023 University of Utah, The Trustees of Columbia University in
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

#include "AkimaSpline.h"
#include <limits>

//--------------------------------------------------------------------------------
struct AkimaSpline::Impl
{
    int                 ncoef;    //! number of B-spline coefficients
    std::vector<double> xknot;    //! knot sequence
    std::vector<double> acoef;    //! Akima spline coefficient a
    std::vector<double> bcoef;    //! Akima spline coefficient b
    std::vector<double> ccoef;    //! Akima spline coefficient c
    std::vector<double> dcoef;    //! Akima spline coefficient d
};

//--------------------------------------------------------------------------------
AkimaSpline::AkimaSpline() : im(new AkimaSpline::Impl)
{
    im->ncoef = 0;
}

//--------------------------------------------------------------------------------
// destructor
AkimaSpline::~AkimaSpline()
{
    delete im;
    im = nullptr;
}

//--------------------------------------------------------------------------------
// copy constructor
AkimaSpline::AkimaSpline(const AkimaSpline& as) : im(new AkimaSpline::Impl)
{
    im->ncoef  = as.im->ncoef;
    im->xknot  = as.im->xknot;
    im->acoef  = as.im->acoef;
    im->bcoef  = as.im->bcoef;
    im->ccoef  = as.im->ccoef;
    im->dcoef  = as.im->dcoef;
}

//--------------------------------------------------------------------------------
void AkimaSpline::operator = (const AkimaSpline& as)
{
    im->ncoef  = as.im->ncoef;
    im->xknot  = as.im->xknot;
    im->acoef  = as.im->acoef;
    im->bcoef  = as.im->bcoef;
    im->ccoef  = as.im->ccoef;
    im->dcoef  = as.im->dcoef;
}

//--------------------------------------------------------------------------------
// initialize B-spline, using p as control points
bool AkimaSpline::init(const std::vector<vec2d>& p)
{
    int ncoef = (int)p.size();
	if (ncoef < 2) return false;

    im->ncoef = ncoef;
    im->xknot.resize(ncoef);
    im->acoef.resize(ncoef-1);
    im->bcoef.resize(ncoef-1);
    im->ccoef.resize(ncoef-1);
    im->dcoef.resize(ncoef-1);

    // extract m sequence
    std::vector<double> m(ncoef-1);
    for (int i=0; i< ncoef-1; ++i)
        m[i] = (p[i+1].y()-p[i].y())/(p[i+1].x()-p[i].x());
    
    // extract slopes s
    const double eps = 10*std::numeric_limits<double>::epsilon();
    std::vector<double> s(ncoef);
    s[0] = m[0];
    s[1] = (m[0]+m[1])/2;
    if (ncoef > 2) {
        s[ncoef-2] = (m[ncoef-3]+m[ncoef-2])/2;
        s[ncoef-1] = m[ncoef-2];
        for (int i=2; i<ncoef-2; ++i) {
            double d = fabs(m[i+1]-m[i]) + fabs(m[i-1]-m[i-2]);
            s[i] = (fabs(d) > eps) ? (fabs(m[i+1]-m[i])*m[i-1]+fabs(m[i-1]-m[i-2])*m[i])/d : (m[i-1]+m[2])/2;
        }
    }

    // evaluate knots and coefficients
    if (ncoef == 2) {
        double dx = p[1].x()-p[0].x();
        if (fabs(dx) <= eps) return false;
        im->xknot[0] = p[0].x();
        im->xknot[1] = p[1].x();
        im->acoef[0] = p[0].y();
        im->bcoef[0] = s[0];
        im->ccoef[0] = 0;
        im->dcoef[0] = 0;
    }
    else {
        for (int i=0; i<ncoef-1; ++i) {
            double dx = p[i+1].x()-p[i].x();
            if (fabs(dx) <= eps) return false;
            im->xknot[i] = p[i].x();
            im->acoef[i] = p[i].y();
            im->bcoef[i] = s[i];
            im->ccoef[i] = (3*m[i] - 2*s[i] - s[i+1])/dx;
            im->dcoef[i] = (s[i] + s[i+1] - 2*m[i])/(dx*dx);
        }
        im->xknot[ncoef-1] = p[ncoef-1].x();
    }
    
    return true;
}

//--------------------------------------------------------------------------------
// evaluate Akima spline at x
double AkimaSpline::eval(double x) const
{
    // perform binary search to locate knot interval that encloses x
    int j = 0, jh = im->ncoef-1;
    while (jh - j > 1) {
        int jm = (j+jh)/2;
        if ((im->xknot[j] <= x) && (x < im->xknot[jm]))
            jh = jm;
        else
            j = jm;
    }
    
    double dx = x - im->xknot[j];
    double y = im->acoef[j] + dx*(im->bcoef[j] + dx*(im->ccoef[j] + dx*im->dcoef[j]));
    
    return y;
}

//--------------------------------------------------------------------------------
double AkimaSpline::eval_deriv(double x) const
{
    // perform binary search to locate knot interval that encloses x
    int j = 0, jh = im->ncoef-1;
    while (jh - j > 1) {
        int jm = (j+jh)/2;
        if ((im->xknot[j] <= x) && (x < im->xknot[jm]))
            jh = jm;
        else
            j = jm;
    }
    
    double dx = x - im->xknot[j];
    double dy = im->bcoef[j] + dx*(2*im->ccoef[j] + 3*dx*im->dcoef[j]);
    
    return dy;
}

//--------------------------------------------------------------------------------
double AkimaSpline::eval_deriv2(double x) const
{
    // perform binary search to locate knot interval that encloses x
    int j = 0, jh = im->ncoef-1;
    while (jh - j > 1) {
        int jm = (j+jh)/2;
        if ((im->xknot[j] <= x) && (x < im->xknot[jm]))
            jh = jm;
        else
            j = jm;
    }
    
    double dx = x - im->xknot[j];
    double d2y = 2*(im->ccoef[j] + 3*dx*im->dcoef[j]);
    
    return d2y;
}
