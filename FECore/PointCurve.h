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
#include "fecore_api.h"
#include "vec2d.h"
#include <vector>

class FECORE_API PointCurve
{
	class Imp;

public:
	//! Interpolation functions
	enum INTFUNC { LINEAR = 0, STEP = 1, SMOOTH = 2, CSPLINE = 3, CPOINTS = 4, APPROX = 5, SMOOTH_STEP = 6 };

	//! Extend mode
	enum EXTMODE { CONSTANT, EXTRAPOLATE, REPEAT, REPEAT_OFFSET };

public:
	//! default constructor
	PointCurve();

	//! copy constructor
	PointCurve(const PointCurve& pc);

	//! assignment operator
	void operator = (const PointCurve& pc);

	//! destructor
	~PointCurve();

	//! call this to update internal data structures
	bool Update();

	//! adds a point to the point curve
	int Add(double x, double y);

	//! adds a point to the point curve
	int Add(const vec2d& p);

	//! Clears the loadcurve data
	void Clear();

	//! set the x and y value of point i
	void SetPoint(int i, double x, double y);
	void SetPoint(int i, const vec2d& p);

	//! set all points at once
	void SetPoints(const std::vector<vec2d>& points);

	//! return all points
	std::vector<vec2d> GetPoints() const;

	//! remove a point
	void Delete(int n);

	//! remove several points at once
	void Delete(const std::vector<int>& indexList);

	//! Set the type of interpolation
	void SetInterpolator(int fnc);

	//! return current interpolator
	int GetInterpolator() const;

	//! Set the extend mode
	void SetExtendMode(int mode);

	//! Get the extend mode
	int GetExtendMode() const;

	//! get a point
	vec2d Point(int i) const;

	//! finds closest load point
	int FindPoint(double t, double& tval, int startIndex = 0);

	//! return nr of points
	int Points() const;

	//! see if there is a point at time t
	bool HasPoint(double t) const;

public: // operations
	
	// scale all y points by s
	void Scale(double s);

public:

	//! returns the value of the load curve at time
	double value(double x) const;

	//! returns the derivative value at time
	double derive(double x) const;

	//! returns the second derivative value at time
	double deriv2(double x) const;

	//! returns the definite integral value between a and b
	double integrate(double a, double b) const;

protected:
	double ExtendValue(double t) const;

private:
	Imp* im;
};
