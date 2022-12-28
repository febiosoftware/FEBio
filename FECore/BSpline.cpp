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

#include "BSpline.h"
#include "matrix.h"
#include <limits>

//--------------------------------------------------------------------------------
struct BSpline::Impl 
{
    int                 korder;   //! B-spline order
    int                 ncoef;    //! number of B-spline coefficients
    std::vector<double> xknot;    //! knot sequence
    std::vector<double> coeff;    //! B-spline coefficients
};

//--------------------------------------------------------------------------------
BSpline::BSpline() : im(new BSpline::Impl)
{
    im->korder = 0;
    im->ncoef = 0;
}

//--------------------------------------------------------------------------------
// destructor
BSpline::~BSpline()
{
    delete im;
    im = nullptr;
}

//--------------------------------------------------------------------------------
// copy constructor
BSpline::BSpline(const BSpline& bs) : im(new BSpline::Impl)
{
    im->korder = bs.im->korder;
    im->ncoef  = bs.im->ncoef;
    im->xknot  = bs.im->xknot;
    im->coeff  = bs.im->coeff;
}

//--------------------------------------------------------------------------------
void BSpline::operator = (const BSpline& bs)
{
    im->korder = bs.im->korder;
    im->ncoef  = bs.im->ncoef;
    im->xknot  = bs.im->xknot;
    im->coeff  = bs.im->coeff;
}

//--------------------------------------------------------------------------------
// initialize B-spline, using p as control points
bool BSpline::init(int korder, const std::vector<vec2d>& p)
{
    int ncoef = (int)p.size();
	if (ncoef < 2) return false;

    im->korder = korder;
    im->ncoef = ncoef;
    if (ncoef < korder) return false;
    im->coeff.resize(ncoef);
    
    // extract breakpoint sequence and spline coefficients
    std::vector<double> bp(ncoef);
    for (int i=0; i< ncoef; ++i) {
        bp[i] = p[i].x();
        im->coeff[i] = p[i].y();
    }
    
    // evaluate knot sequence from breakpoint sequence
    im->xknot.resize(korder + ncoef);
    int khalf;
    const double eps = 10*std::numeric_limits<double>::epsilon();
    double e = eps*(bp[ncoef-1] - bp[ncoef-2]);
    
    if ((korder % 2) == 0) {
        khalf = korder/2;
        for (int i=0; i< korder; ++i)
            im->xknot[i] = bp.front();
        for (int i= korder; i< ncoef; ++i)
            im->xknot[i] = bp[i-khalf];
        for (int i=ncoef; i<ncoef+korder; ++i)
            im->xknot[i]  = bp.back() + e;
    }
    else {
        khalf = (korder-1)/2;
        for (int i=0; i<korder; ++i)
            im->xknot[i] = bp.front();
        for (int i=korder; i<ncoef; ++i)
            im->xknot[i] = (bp[i-khalf] + bp[i-1-khalf])/2;
        for (int i=ncoef; i<ncoef+korder; ++i)
            im->xknot[i]  = bp.back() + e;
    }

    return true;
}

//--------------------------------------------------------------------------------
// evaluate B-spline at x using de Boor algorithm (de Boor 1986)
double BSpline::eval(double x, int korder, const std::vector<double>& xknot,
            int ncoef, const std::vector<double>& coeff) const
{
    // perform binary search to locate knot interval that encloses x
    int j = korder-1, jh = ncoef;
    while (jh - j > 1) {
        int jm = (j+jh)/2;
        if ((xknot[j] <= x) && (x < xknot[jm]))
            jh = jm;
        else
            j = jm;
    }
    
    std::vector<double> c = coeff;
    double w;
    for (int r=0; r<korder-1; ++r) {
        for (int i=j; i>j-korder+r+1; --i) {
            if (xknot[i] != xknot[i+korder-r-1])
                w = (x - xknot[i])/(xknot[i+korder-r-1] - xknot[i]);
            else
                w = 0;
            c[i] = (1-w)*c[i-1] + w*c[i];
        }
    }
    return c[j];
}

//--------------------------------------------------------------------------------
// evaluate B-spline at x using de Boor algorithm (de Boor 1986)
double BSpline::eval(double x) const
{
    return eval(x, im->korder, im->xknot, im->ncoef, im->coeff);
}

//--------------------------------------------------------------------------------
double BSpline::eval_deriv(double x) const { return eval_nderiv(x, 1); }

//--------------------------------------------------------------------------------
double BSpline::eval_deriv2(double x) const { return eval_nderiv(x, 2); }

//--------------------------------------------------------------------------------
// evaluate B-spline n-th derivative at x using de Boor algorithm (de Boor 1986)
double BSpline::eval_nderiv(double x, int n) const
{
    double deriv = 0;

    int korder = im->korder;
    int ncoef = im->ncoef;
    if (n >= korder) return deriv;
    
    std::vector<double> coeff = im->coeff;

    for (int k=1; k<=n; ++k) {
        for (int i= im->ncoef-1; i>=k; --i) {
            if (im->xknot[i+korder-k] > im->xknot[i]) {
                coeff[i] = (korder - k)*(coeff[i] - coeff[i-1])/
                (im->xknot[i+korder-k] - im->xknot[i]);
            }
        }
    }

    // extract portion of vectors
    std::vector<double> xknot(im->xknot.begin() + n, im->xknot.end());
    std::vector<double> doeff(coeff.begin() + n,coeff.end());

    deriv = eval(x,korder-n, xknot, ncoef-n, doeff);
    
    return deriv;
}

//--------------------------------------------------------------------------------
// evaluate B-spline blending functions at x
std::vector<double> BSpline::blending_functions(double x) const
{
    int korder = im->korder;
    int ncoef  = im->ncoef;

    std::vector<double> bsbldg(ncoef);
    std::vector<std::vector<double> > d(ncoef, std::vector<double>(korder));
    
    for (int i=0; i<ncoef; ++i) {
        if ((im->xknot[i] <= x) && (x < im->xknot[i+1]))
            d[i][0] = 1;
        else
            d[i][0] = 0;
    }

    double r1, r2;
    for (int k=2; k<=korder; ++k) {
        for (int i=0; i<ncoef; ++i) {
            if (im->xknot[i+k-1] == im->xknot[i])
                r1 = 0;
            else
                r1 = (x- im->xknot[i])*d[i][k-2]/(im->xknot[i+k-1]-im->xknot[i]);
            if (im->xknot[i+k] == im->xknot[i+1])
                r2 = 0;
            else
                r2 = (im->xknot[i+k]-x)*d[i+1][k-2]/(im->xknot[i+k]- im->xknot[i+1]);
            d[i][k-1] = r1 + r2;
        }
    }
    
    for (int i=0; i<ncoef; ++i)
        bsbldg[i] = d[i][korder-1];
    
    return bsbldg;
}

//--------------------------------------------------------------------------------
// fit a B-spline of order korder, with ncoef coefficients, to the points p
bool BSpline::fit(int korder, int ncoef, const std::vector<vec2d>& p)
{
    // check for valid spline order
    if (korder <1) return false;
    
    // number of points to fit
    int np = (int)p.size();
    if (np < korder) return false;
    
    // for an interpolation, use p.x as breakpoints to generate knot sequence
    bool binit = true;
    if (ncoef == np) binit = init(korder, p);
    // otherwise, generate breakpoints uniformly over range of x
    else {
        std::vector<vec2d> q(ncoef,vec2d(0, 0));
        double dx = (p[np-1].x() - p[0].x())/(ncoef-1);
        for (int i=0; i<ncoef; ++i)
            q[i].x() = p[0].x() + i*dx;
        binit = init(korder, q);
    }
    if (binit == false) return binit;
    
    // evaluate B-spline blending functions at p.x using this knot sequence
    matrix wk1(np,ncoef);
    for (int j=0; j<np; ++j) {
        std::vector<double> bsbldg = blending_functions(p[j].x());
        for (int i=0; i<ncoef; ++i) wk1(j,i) = bsbldg[i];
    }
    
    // evaluate the coefficient matrix
    matrix wk2 = wk1.transpose()*wk1;
    
    // evaluate the right-hand-side
    std::vector<double> rhs(ncoef,0);
    for (int k=0; k<ncoef; ++k)
        for (int j=0; j<np; ++j)
            rhs[k] += p[j].y()*wk1(j,k);

    // solve the system of equations
    wk2.solve(im->coeff, rhs);
    
    return true;
}

//--------------------------------------------------------------------------------
// use given points as interpolation points
bool BSpline::init_interpolation(int korder, const std::vector<vec2d>& p)
{
    int ncoef = (int)p.size();
    return fit(korder, ncoef, p);
}

//--------------------------------------------------------------------------------
// perform spline approximation over points p, using ncoef coefficients
bool BSpline::init_approximation(int korder, int ncoef, const std::vector<vec2d>& p)
{
    return fit(korder, ncoef, p);
}
