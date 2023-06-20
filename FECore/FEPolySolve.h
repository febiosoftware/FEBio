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

#pragma once
#include "fecore_api.h"
#include <complex>
#include <vector>

class FECORE_API FEPolySolve
{
public:
    // constructor
    FEPolySolve() {}
    
    void fn(std::complex<double>& z, std::complex<double>& fz, std::vector<double> a);
    
    bool dflate(std::complex<double> zero, const int i, int& kount,
                std::complex<double>& fzero, std::complex<double>& fzrdfl,
                std::complex<double>* zeros, std::vector<double> a);
    
    bool muller(bool fnreal, std::complex<double>* zeros, const int n, const int nprev,
                const int maxit, const double ep1, const double ep2, std::vector<double> a);
    
    bool newton(double& zero, const int n, const int maxit,
                const double ep1, const double ep2, std::vector<double> a);
    
    bool poly1(std::vector<double> a, double& x);

    bool poly2(std::vector<double> a, double& x);
    
    bool polyn(int n, std::vector<double> a, double& x);
    
    bool polym(int n, std::vector<double> a, double& x);
    
    bool solvepoly(int n, std::vector<double> a, double& x, bool nwt=true);
};
