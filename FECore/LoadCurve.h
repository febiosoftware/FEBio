#pragma once
#include "FECoreBase.h"
#include "FEFunction1D.h"

//-----------------------------------------------------------------------------
// Base class for load curves.
// Load curves are used to manipulate the time dependency of model parameters.
class FECORE_API FELoadCurve
{
public:
	// constructor
	FELoadCurve();
	FELoadCurve(FEFunction1D* fnc);
	FELoadCurve(const FELoadCurve& lc);

	void operator = (const FELoadCurve& lc);

	// assign a function
	void SetFunction(FEFunction1D* f);

	// destructor
	virtual ~FELoadCurve();

	// return the last evaluated function value
	double Value() const { return m_value; }

	//! evaluates the loadcurve at time
	void Evaluate(double time)
	{
		m_value = m_fnc->value(time);
	}

	// evaluate a load curve
	double Value(double time)
	{
		return m_fnc->value(time);
	}

	void Serialize(DumpStream& ar);

	bool CopyFrom(FELoadCurve* lc);

	// evaluate the derivative at time t
	double Deriv(double t) const { return m_fnc->derive(t); }

	FEFunction1D* GetFunction() { return m_fnc; }

private:
	double	m_value;	//!< value of last call to Value

	FEFunction1D*	m_fnc;	//!< functin to evaluate
};
