#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
// Base class for load curves.
class FELoadCurve : public FECoreBase
{
public:
	// constructor
	FELoadCurve();
	FELoadCurve(const FELoadCurve& lc);

	void operator = (const FELoadCurve& lc);

	// destructor
	virtual ~FELoadCurve();

	// return the last evaluated function value
	double Value() const { return m_value; }

	//! evaluates the loadcurve at time
	void Evaluate(double time)
	{
		m_value = Value(time);
	}

	void Serialize(DumpStream& ar);

public:
	// evaluate the function at time t
	virtual double Value(double t) const = 0;

	// evaluate the derivative at time t
	virtual double Deriv(double t) const = 0;

private:
	double	m_value;	//!< value of last call to Value
};

//-----------------------------------------------------------------------------
// A loadcurve that generates a linear ramp
class FELinearRamp : public FELoadCurve
{
public:
	FELinearRamp(FEModel* fem) : m_slope(0.0), m_intercept(0.0) {}
	FELinearRamp(double m, double y0) : m_slope(m), m_intercept(y0){}

	double Value(double t) const
	{
		return m_slope*t  + m_intercept;
	}

	double Deriv(double t) const
	{
		return m_slope;
	}

private:
	double	m_slope;
	double	m_intercept;
};
