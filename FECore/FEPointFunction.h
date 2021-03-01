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
#include "FEFunction1D.h"
#include <vector>

//-----------------------------------------------------------------------------
class DumpStream;

//-----------------------------------------------------------------------------
//! This class implements a function defined by a set of ordered (point,value) pairs.
//! It uses an interpolation scheme to interpolate between data points and also 
//! has a way of specifying the function value outside of the domain defined by the first and 
//! last data point. 

class FECORE_API FEPointFunction : public FEFunction1D
{
	class Imp;

public:
	//! Load point structure
	struct LOADPOINT
	{
		double time;
		double value;
	};

public:
	//! Interpolation functions
	enum INTFUNC { STEP = 0, LINEAR = 1, SMOOTH = 2, POLYNOMIAL = 3, CSPLINE = 4, AKIMA = 5, STEFFEN = 6 };

	//! Extend mode
	enum EXTMODE { CONSTANT, EXTRAPOLATE, REPEAT, REPEAT_OFFSET };

public:
	//! default constructor
	FEPointFunction(FEModel* fem);

	//! destructor
	~FEPointFunction();
    
    //! initialize
    bool Init() override;

	//! adds a point to the point curve
	void Add(double x, double y);

	//! Clears the loadcurve data
	void Clear() override;

	//! set the x and y value of point i
	void SetPoint(int i, double x, double y);

	//! Set the type of interpolation
	void SetInterpolation(INTFUNC fnc) { m_fnc = fnc; }

	//! Set the extend mode
	void SetExtendMode(EXTMODE mode) { m_ext = mode; }

	//! returns point i
	LOADPOINT LoadPoint(int i) const;

	//! finds closest load point
	int FindPoint(double t, double& tval, int startIndex = 0);

	//! return nr of points
	int Points() const;

	//! see if there is a point at time t
	bool HasPoint(double t) const;

	//! Serialize data to archive
	void Serialize(DumpStream& ar) override;

	// copy data from other curve
	FEFunction1D* copy() override;

	// copy from another function
	void CopyFrom(const FEPointFunction& f);

public: // operations

	// scale all y points by s
	void Scale(double s);

public: // implement from base class

	//! returns the value of the load curve at time
	double value(double x) const override;

	//! returns the derivative value at time
	double derive(double x) const override;

    //! returns the second derivative value at time
    double deriv2(double x) const override;

	//! returns the definite integral value between a and b
	double integrate(double a, double b) const override;

protected:
	double ExtendValue(double t) const;

// private:
// 	//! returns the area of a trapezoid between a and b
// 	double trap(double a, double b) const;

	// TODO: I need to make this public so the parameters can be mapped to the FELoadCurve
public:
	int		m_fnc;	//!< interpolation function
	int		m_ext;	//!< extend mode
    bool    m_bln;  //!< points represent (ln(x),y) instead of (x,y)
	std::vector<vec2d>	m_points;
    
private:
	Imp*	imp;

	DECLARE_FECORE_CLASS();
};

typedef FEPointFunction::LOADPOINT LOADPOINT;
