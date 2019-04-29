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
#include "fecore_api.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEModel;
class DumpStream;

//-----------------------------------------------------------------------------
// Class that represents a 1D function
// Derived classes need to implement:
//   - value(double): This evaluates the function at a particular point
//   - copy()       : Create a copy of the class
//   - derive(double) : Calculate the derivative. This is optional, but allows implementation of more efficient algorithm, since default implements forward difference
//
class FECORE_API FEFunction1D : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FEFunction1D(FEModel* pfem);

	// serialization
	void Serialize(DumpStream& ar);

	// this class requires a copy member
	virtual FEFunction1D* copy() = 0;

	// evaluate the function at x
	// must be defined by derived classes
	virtual double value(double x) const = 0;

	// value of first derivative of function at x
	// can be overridden by derived classes.
	// default implementation is a forward-difference
	virtual double derive(double x) const;

	virtual void Clear() {}
};

//-----------------------------------------------------------------------------
// A linear function
class FECORE_API FELinearFunction : public FEFunction1D
{
public:
	FELinearFunction(FEModel* fem) : FEFunction1D(fem), m_slope(0.0), m_intercept(0.0) {}
	FELinearFunction(FEModel* fem, double m, double y0) : FEFunction1D(fem), m_slope(m), m_intercept(y0) {}
	FEFunction1D* copy() { return new FELinearFunction(GetFEModel(), m_slope, m_intercept); }

	double value(double t) const
	{
		return m_slope*t + m_intercept;
	}

	double deriv(double t) const
	{
		return m_slope;
	}

private:
	double	m_slope;
	double	m_intercept;
};
