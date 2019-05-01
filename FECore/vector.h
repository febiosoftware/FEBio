/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include <math.h>
#include <memory.h>
#include <vector>
#include <algorithm>
#include "vec3d.h"
#include "fecore_api.h"
using namespace std;

class FEMesh;
class FEDofList;

double FECORE_API operator*(const vector<double>& a, const vector<double>& b);
vector<double> FECORE_API operator - (vector<double>& a, vector<double>& b);
template<typename T> void zero(vector<T>& a) { fill(a.begin(), a.end(), T(0)); }
template<> inline void zero<vec3d>(vector<vec3d>& a) { fill(a.begin(), a.end(), vec3d(0,0,0)); }
template<typename T> void assign(vector<T>& a, const T& v) { fill(a.begin(), a.end(), v); }
void FECORE_API operator+=(vector<double>& a, const vector<double>& b);
void FECORE_API operator-=(vector<double>& a, const vector<double>& b);
void FECORE_API operator*=(vector<double>& a, double b);
vector<double> FECORE_API operator+(const vector<double>& a, const vector<double>& b);

// copy vector and scale
void FECORE_API vcopys(vector<double>& a, const vector<double>& b, double s);

// add scaled vector
void FECORE_API vadds(vector<double>& a, const vector<double>& b, double s);
void FECORE_API vsubs(vector<double>& a, const vector<double>& b, double s);

// vector subtraction: a = l - r
void FECORE_API vsub(vector<double>& a, const vector<double>& l, const vector<double>& r);

// scale each component of a vector
void FECORE_API vscale(vector<double>& a, const vector<double>& s);

// gather operation (copy mesh data to vector)
void FECORE_API gather(vector<double>& v, FEMesh& mesh, int ndof);
void FECORE_API gather(vector<double>& v, FEMesh& mesh, const vector<int>& dof);

// scatter operation (copy vector data to mesh)
void FECORE_API scatter(vector<double>& v, FEMesh& mesh, int ndof);
void FECORE_API scatter3(vector<double>& v, FEMesh& mesh, int ndof1, int ndof2, int ndof3);
void FECORE_API scatter(vector<double>& v, FEMesh& mesh, const FEDofList& dofs);

// calculate l2 norm of vector
double FECORE_API l2_norm(const vector<double>& v);
double l2_norm(double* x, int n);
