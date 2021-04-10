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


#include "BSpline.h"
#include "matrix.h"
#include <limits>

//--------------------------------------------------------------------------------
BSpline::BSpline()
{
    m_korder = 0;
    m_ncoef = 0;
}

//--------------------------------------------------------------------------------
// copy constructor
BSpline::BSpline(const BSpline& bs)
{
    m_korder = bs.m_korder;
    m_ncoef = bs.m_ncoef;
    m_xknot = bs.m_xknot;
    m_coeff = bs.m_coeff;
}

//--------------------------------------------------------------------------------
// initialize B-spline, using p as control points
bool BSpline::init(int korder, std::vector<vec2d> p)
{
    m_korder = korder;
    m_ncoef = (int)p.size();
    if (m_ncoef < m_korder) return false;
    m_coeff.resize(m_ncoef);
    
    // extract breakpoint sequence and spline coefficients
    std::vector<double> bp(m_ncoef);
    for (int i=0; i<m_ncoef; ++i) {
        bp[i] = p[i].x();
        m_coeff[i] = p[i].y();
    }
    
    // evaluate knot sequence from breakpoint sequence
    m_xknot.resize(m_korder + m_ncoef);
    int khalf;
    const double eps = 10*std::numeric_limits<double>::epsilon();
    double e = eps*(bp[m_ncoef-1] - bp[m_ncoef-2]);
    
    if ((m_korder % 2) == 0) {
        khalf = m_korder/2;
        for (int i=0; i<m_korder; ++i)
            m_xknot[i] = bp.front();
        for (int i=m_korder; i<m_ncoef; ++i)
            m_xknot[i] = bp[i-khalf];
        for (int i=m_ncoef; i<m_ncoef+m_korder; ++i)
            m_xknot[i]  = bp.back() + e;
    }
    else {
        khalf = (m_korder-1)/2;
        for (int i=0; i<m_korder; ++i)
            m_xknot[i] = bp.front();
        for (int i=m_korder; i<m_ncoef; ++i)
            m_xknot[i] = (bp[i-khalf] + bp[i-1-khalf])/2;
        for (int i=m_ncoef; i<m_ncoef+m_korder; ++i)
            m_xknot[i]  = bp.back() + e;
    }

    return true;
}

//--------------------------------------------------------------------------------
// evaluate B-spline at x using de Boor algorithm (de Boor 1986)
double BSpline::eval(double x, int korder, std::vector<double>xknot,
            int ncoef, std::vector<double>coeff)
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
    
    double w;
    for (int r=0; r<korder-1; ++r) {
        for (int i=j; i>j-korder+r+1; --i) {
            if (xknot[i] != xknot[i+korder-r-1])
                w = (x - xknot[i])/(xknot[i+korder-r-1] - xknot[i]);
            else
                w = 0;
            coeff[i] = (1-w)*coeff[i-1] + w*coeff[i];
        }
    }
    return coeff[j];
}

//--------------------------------------------------------------------------------
// evaluate B-spline at x using de Boor algorithm (de Boor 1986)
double BSpline::eval(double x)
{
    std::vector<double> coeff = m_coeff;
    return eval(x,m_korder, m_xknot, m_ncoef, coeff);
}

//--------------------------------------------------------------------------------
// evaluate B-spline n-th derivative at x using de Boor algorithm (de Boor 1986)
double BSpline::eval_nderiv(double x, int n)
{
    double deriv = 0;
    
    if (n >= m_korder) return deriv;
    
    std::vector<double> coeff = m_coeff;

    for (int k=1; k<=n; ++k) {
        for (int i=m_ncoef-1; i>=k; --i) {
            if (m_xknot[i+m_korder-k] > m_xknot[i]) {
                coeff[i] = (m_korder - k)*(coeff[i] - coeff[i-1])/
                (m_xknot[i+m_korder-k] - m_xknot[i]);
            }
        }
    }

    // extract portion of vectors
    std::vector<double> xknot(m_xknot.begin() + n,m_xknot.end());
    std::vector<double> doeff(coeff.begin() + n,coeff.end());

    deriv = eval(x,m_korder-n, xknot, m_ncoef-n, doeff);
    
    return deriv;
}

//--------------------------------------------------------------------------------
// evaluate B-spline blending functions at x
std::vector<double> BSpline::blending_functions(double x)
{
    std::vector<double> bsbldg(m_ncoef);
    std::vector<std::vector<double> > d(m_ncoef, std::vector<double>(m_korder));
    
    for (int i=0; i<m_ncoef; ++i) {
        if ((m_xknot[i] <= x) && (x < m_xknot[i+1]))
            d[i][0] = 1;
        else
            d[i][0] = 0;
    }

    double r1, r2;
    for (int k=2; k<=m_korder; ++k) {
        for (int i=0; i<m_ncoef; ++i) {
            if (m_xknot[i+k-1] == m_xknot[i])
                r1 = 0;
            else
                r1 = (x-m_xknot[i])*d[i][k-2]/(m_xknot[i+k-1]-m_xknot[i]);
            if (m_xknot[i+k] == m_xknot[i+1])
                r2 = 0;
            else
                r2 = (m_xknot[i+k]-x)*d[i+1][k-2]/(m_xknot[i+k]-m_xknot[i+1]);
            d[i][k-1] = r1 + r2;
        }
    }
    
    for (int i=0; i<m_ncoef; ++i)
        bsbldg[i] = d[i][m_korder-1];
    
    return bsbldg;
}

//--------------------------------------------------------------------------------
// fit a B-spline of order korder, with ncoef coefficients, to the points p
bool BSpline::fit(int korder, int ncoef, std::vector<vec2d>p)
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
    wk2.solve(rhs, m_coeff);
    
    return true;
}

//--------------------------------------------------------------------------------
// use given points as interpolation points
bool BSpline::init_interpolation(int korder, std::vector<vec2d> p)
{
    int ncoef = (int)p.size();
    return fit(korder, ncoef, p);
}

//--------------------------------------------------------------------------------
// perform spline approximation over points p, using ncoef coefficients
bool BSpline::init_approximation(int korder, int ncoef, std::vector<vec2d> p)
{
    return fit(korder, ncoef, p);
}
