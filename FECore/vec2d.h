/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "fecore_api.h"

class vec2d
{
public:
	// constructor
	vec2d() { r[0] = r[1] = 0; }
	explicit vec2d(double v) { r[0] = r[1] = v; }
	vec2d(double x, double y) { r[0] = x; r[1] = y; }

	bool operator == (const vec2d& r) const { return (r[0] == r.r[0]) && (r[1] == r.r[1]); }

	// access operators
	double operator [] (int i) const { return r[i]; }
	double& operator [] (int i) { return r[i]; }

	double& x() { return r[0]; }
	double& y() { return r[1]; }

	double x() const { return r[0]; }
	double y() const { return r[1]; }

	double norm() const { return sqrt(r[0] * r[0] + r[1] * r[1]); }
	double norm2() const { return r[0] * r[0] + r[1] * r[1]; }

	double unit() { double R = norm(); if (R != 0) { r[0] /= R; r[1] /= R; }; return R; }

public: // arithmetic operators

	vec2d operator + (const vec2d& v) const { return vec2d(r[0]+v.r[0], r[1]+v.r[1]); }
	vec2d operator - (const vec2d& v) const { return vec2d(r[0]-v.r[0], r[1]-v.r[1]); }
	vec2d operator * (double g) const { return vec2d(r[0]*g, r[1]*g); }
	vec2d operator / (double g) const { return vec2d(r[0]/g, r[1]/g); }

	vec2d& operator += (const vec2d& v) { r[0] += v.r[0]; r[1] += v.r[1]; return *this; }
	vec2d& operator -= (const vec2d& v) { r[0] -= v.r[0]; r[1] -= v.r[1]; return *this; }
	vec2d& operator *= (double g) { r[0] *= g; r[1] *= g; return *this; }
	vec2d& operator /= (double g) { r[0] /= g; r[1] /= g; return *this; }

    vec2d operator - () const { return vec2d(-r[0], -r[1]); }
    
	// dot product
	double operator * (const vec2d& v) const { return r[0]*v[0] + r[1]*v[1]; }

public:
	double r[2];
};

//-----------------------------------------------------------------------------
class vec2i
{
public:
	vec2i() { x = y = 0; }
	vec2i(int X, int Y) { x = X; y = Y; }

    bool operator == (const vec2i& r) const { return (x == r.x) && (y == r.y); }

public:
	int		x, y;
};

//-----------------------------------------------------------------------------
class vec2f
{
public:
	vec2f() { x = y = 0.f; }
	vec2f(float rx, float ry) { x = rx; y = ry; }

public:
	float	x, y;
};
