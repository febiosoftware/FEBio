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

#include "FEPolySolve.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
//! Polynomial root solver

// function whose roots needs to be evaluated
void FEPolySolve::fn(std::complex<double>& z, std::complex<double>& fz, std::vector<double> a)
{
    int n = (int)a.size()-1;
    fz = a[0];
    std::complex<double> x(1,0);
    
    for (int i=1; i<=n; ++i) {
        x *= z;
        fz += a[i]*x;
    }
    return;
}

//-----------------------------------------------------------------------------
// deflation
bool FEPolySolve::dflate(std::complex<double> zero, const int i, int& kount,
                         std::complex<double>& fzero, std::complex<double>& fzrdfl,
                         std::complex<double>* zeros, std::vector<double> a)
{
    std::complex<double> den;
    ++kount;
    fn(zero, fzero, a);
    fzrdfl = fzero;
    if (i < 1) return false;
    for (int j=0; j<i; ++j) {
        den = zero - zeros[j];
        if (abs(den) == 0) {
            zeros[i] = zero*1.001;
            return true;
        } else {
            fzrdfl = fzrdfl/den;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
// Muller's method for solving roots of a function
bool FEPolySolve::muller(bool fnreal, std::complex<double>* zeros, const int n, const int nprev,
            const int maxit, const double ep1, const double ep2, std::vector<double> a)
{
    int kount;
    std::complex<double> dvdf1p, fzrprv, fzrdfl, divdf1, divdf2;
    std::complex<double> fzr, zero, c, den, sqr, z;
    
    // initialization
    double eps1 = (ep1 > 1e-12) ? ep1:1e-12;
    double eps2 = (ep2 > 1e-20) ? ep2:1e-20;
    
    for (int i=nprev; i<n; ++i) {
        kount = 0;
    eloop:
        zero = zeros[i];
        std::complex<double> h = 0.5;
        std::complex<double> hprev = -1.0;
        
        // compute first three estimates for zero as
        // zero+0.5, zero-0.5, zero
        z = zero + 0.5;
        if (dflate(z, i, kount, fzr, dvdf1p, zeros, a)) goto eloop;
        z = zero - 0.5;
        if (dflate(z, i, kount, fzr, fzrprv, zeros, a)) goto eloop;
        dvdf1p = (fzrprv - dvdf1p)/hprev;
        if (dflate(zero, i, kount, fzr, fzrdfl, zeros, a)) goto eloop;
        do {
            divdf1 = (fzrdfl - fzrprv)/h;
            divdf2 = (divdf1 - dvdf1p)/(h+hprev);
            hprev = h;
            dvdf1p = divdf1;
            c = divdf1 + h*divdf2;
            sqr = c*c - 4.*fzrdfl*divdf2;
            if (fnreal && (sqr.real() < 0)) sqr = 0;
            sqr = sqrt(sqr);
            if (c.real()*sqr.real()+c.imag()*sqr.imag() < 0) {
                den = c - sqr;
            } else {
                den = c + sqr;
            }
            if (abs(den) <= 0.) den = 1.;
            h = -2.*fzrdfl/den;
            fzrprv = fzrdfl;
            zero = zero + h;
        dloop:
            fn(zero,fzrdfl,a);
            // check for convergence
            if (abs(h) < eps1*abs(zero)) break;
            if (abs(fzrdfl) < eps2) break;
            // check for divergence
            if (abs(fzrdfl) >= 10.*abs(fzrprv)) {
                h /= 2.;
                zero -= h;
                goto dloop;
            }
        } while (kount < maxit);
        zeros[i] = zero;
    }
    return true;
}

//-----------------------------------------------------------------------------
// Newton's method for finding nearest root of a polynomial
bool FEPolySolve::newton(double& zero, const int n, const int maxit,
            const double ep1, const double ep2, std::vector<double> a)
{
    bool done = false;
    bool conv = false;
    int it = 0;
    double f, df, x, dx, xi;
    x = zero;
    
    while (!done) {
        // Evaluate function and its derivative
        xi = x;
        f = a[0] + a[1]*xi;
        df = a[1];
        for (int i=2; i<=n; ++i) {
            df += i*a[i]*xi;
            xi *= x;
            f += a[i]*xi;
        }
        if (df == 0) break;
        // check absolute convergence and don't update x if met
        if (abs(f) < ep2) {
            done = true;
            conv = true;
            zero = x;
            break;
        }
        // evaluate increment in x
        dx = -f/df;
        x += dx;
        ++it;
        // check relative convergence
        if (abs(dx) < ep1*abs(x)) {
            done = true;
            conv = true;
            zero = x;
        }
        // check iteration count
        else if (it > maxit) {
            done = true;
            zero = x;
        }
    }
    return conv;
}

//-----------------------------------------------------------------------------
// linear
bool FEPolySolve::poly1(std::vector<double> a, double& x)
{
    if (a[1]) {
        x = -a[0]/a[1];
        return true;
    } else {
        return false;
    }
}

//-----------------------------------------------------------------------------
// quadratic
bool FEPolySolve::poly2(std::vector<double> a, double& x)
{
    if (a[2]) {
        x = (-a[1]+sqrt(SQR(a[1])-4*a[0]*a[2]))/(2*a[2]);
        return true;
    } else {
        return poly1(a,x);
    }
}

//-----------------------------------------------------------------------------
// higher order
bool FEPolySolve::polyn(int n, std::vector<double> a, double& x)
{
    //    bool fnreal = true;
    //    vector< complex<double> > zeros(n,complex<double>(1,0));
    int maxit = 100;
    double ep1 = 1e-6;
    double ep2 = 1e-12;
    return newton(x, n, maxit,ep1, ep2, a);
}

//-----------------------------------------------------------------------------
// higher order
bool FEPolySolve::polym(int n, std::vector<double> a, double& x)
{
    bool fnreal = true;
    std::vector< std::complex<double> > zeros(n,std::complex<double>(1,0));
    int maxit = 100;
    double ep1 = 1e-6;
    double ep2 = 1e-12;
    
    muller(fnreal, &zeros[0], n, 0, maxit, ep1, ep2, a);
    for (int i=0; i<n; ++i) {
        if (zeros[i].real() != 0) {
            x = zeros[i].real();
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPolySolve::solvepoly(int n, std::vector<double> a, double& x, bool nwt)
{
    switch (n) {
        case 1:
            return poly1(a, x);
            break;
        case 2:
            return poly2(a, x);
        default:
            if (a[n]) {
                return nwt ? polyn(n, a, x) : polym(n, a, x);
            } else {
                return solvepoly(n-1, a, x);
            }
            break;
    }
}

