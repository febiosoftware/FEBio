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
	DECLARE_SUPER_CLASS(FEFUNCTION1D_ID);

public:
	FEFunction1D(FEModel* pfem);

	// serialization
	void Serialize(DumpStream& ar);

	// return the FEModel
	FEModel* GetFEModel() { return m_fem; }

	// this class requires a copy member
	virtual FEFunction1D* copy() = 0;

	// evaluate the function at x
	// must be defined by derived classes
	virtual double value(double x) const = 0;

	// value of first derivative of function at x
	// can be overridden by derived classes.
	// default implementation is a forward-difference
	virtual double derive(double x) const;

private:
	FEModel*	m_fem;	//!< point to model
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
