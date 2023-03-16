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
#include "FECoreBase.h"
#include "MathObject.h"

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
	FECORE_SUPER_CLASS(FEFUNCTION1D_ID)
	FECORE_BASE_CLASS(FEFunction1D);

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

    virtual double deriv2(double x) const = 0;

	// value of first definite integral of funciton from a to b
	// can be overridden by derived classes.
	// default implementation is trapezoidal rule
	virtual double integrate(double a, double b) const;
    
	virtual void Clear() {}
    
    // invert function
    virtual bool invert(const double f0, double &x);
};

//-----------------------------------------------------------------------------
// A constant function
class FECORE_API FEConstFunction : public FEFunction1D
{
public:
	FEConstFunction(FEModel* fem) : FEFunction1D(fem), m_value(0.0) {}
	FEFunction1D* copy() override { return new FEConstFunction(GetFEModel(), m_value); }

	double value(double t) const override { return m_value;	}
	double derive(double t) const override { return 0.0; }
	double deriv2(double t) const override { return 0.0; }

protected:
	FEConstFunction(FEModel* fem, double val) : FEFunction1D(fem), m_value(val) {}

private:
	double	m_value;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// A linear function
class FECORE_API FELinearFunction : public FEFunction1D
{
public:
	FELinearFunction(FEModel* fem) : FEFunction1D(fem), m_slope(0.0), m_intercept(0.0) {}
	FELinearFunction(FEModel* fem, double m, double y0) : FEFunction1D(fem), m_slope(m), m_intercept(y0) {}
	FEFunction1D* copy() override { return new FELinearFunction(GetFEModel(), m_slope, m_intercept); }

	double value(double t) const override
	{
		return m_slope*t + m_intercept;
	}

	double derive(double t) const override
	{
		return m_slope;
	}

    double deriv2(double t) const override
    {
        return 0;
    }

private:
	double	m_slope;
	double	m_intercept;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// A step function
class FECORE_API FEStepFunction : public FEFunction1D
{
public:
	FEStepFunction(FEModel* fem) : FEFunction1D(fem), m_x0(0.0), m_leftVal(0.0), m_rightVal(1.0) {}
	FEStepFunction(FEModel* fem, double x0, double lv, double rv) : FEFunction1D(fem), m_x0(x0), m_leftVal(lv), m_rightVal(rv) {}
	FEFunction1D* copy() override { return new FEStepFunction(GetFEModel(), m_x0, m_leftVal, m_rightVal); }

	double value(double t) const override
	{
		return (t < m_x0 ? m_leftVal : m_rightVal);
	}

	double derive(double t) const override
	{
		return 0.0;
	}

	double deriv2(double t) const override
	{
		return 0.0;
	}

    // invert function has no unique solution
    bool invert(const double f0, double &x) override { return false; }
    
private:
	double	m_x0;
	double	m_leftVal;
	double	m_rightVal;

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
//! function defined via math expression
class FECORE_API FEMathFunction : public FEFunction1D
{
public:
	FEMathFunction(FEModel* fem);

	bool Init() override;

	void Serialize(DumpStream& ar) override;

	FEFunction1D* copy() override;

	double value(double t) const override;

	double derive(double t) const override;

    double deriv2(double t) const override;

	void SetMathString(const std::string& s);

private:
	void evalParams(std::vector<double>& val, double t) const;

	bool BuildMathExpressions();

private:
	std::string			m_s;
	int					m_ix;			// index of independent variable
	std::vector<FEParamValue>	m_var;	// list of model parameters that are used as variables in expression.

	MSimpleExpression	m_exp;
	MSimpleExpression	m_dexp;
    MSimpleExpression   m_d2exp;

	DECLARE_FECORE_CLASS();
};
