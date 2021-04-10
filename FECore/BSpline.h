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



#pragma once
#include "fecore_api.h"
#include "vec2d.h"
#include <vector>

class FECORE_API BSpline
{
public:
    // constructor
    BSpline();
    // copy constructor
    BSpline(const BSpline& bs);
    // destructor
    ~BSpline() {}

public:
    // spline function evaluations
    double eval(double x, int korder, std::vector<double>xknot,
                int ncoef, std::vector<double>coeff);
    double eval(double x);
    double eval_nderiv(double x, int n);
    double eval_deriv(double x) { return eval_nderiv(x, 1); }
    double eval_deriv2(double x) { return eval_nderiv(x, 2); }
    
public:
    // spline fitting
    std::vector<double> blending_functions(double x);
    bool fit(int korder, int ncoef, std::vector<vec2d>p);
    
public:
    // initializations
    // use given points p as control points
    bool init(int korder, std::vector<vec2d> p);
    // use given points p as interpolation points
    bool init_interpolation(int korder, std::vector<vec2d> p);
    // perform spline approximation over points p, using ncoef coefficients
    bool init_approximation(int korder, int ncoef, std::vector<vec2d> p);

private:
    int                 m_korder;   //! B-spline order
    int                 m_ncoef;    //! number of B-spline coefficients
    std::vector<double> m_xknot;    //! knot sequence
    std::vector<double> m_coeff;    //! B-spline coefficients
    
};
