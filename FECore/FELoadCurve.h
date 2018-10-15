#pragma once
#include "FELoadController.h"
#include "FEPointFunction.h"

//-----------------------------------------------------------------------------
// Base class for load curves.
// Load curves are used to manipulate the time dependency of model parameters.
class FECORE_API FELoadCurve : public FELoadController
{
public:
	// constructor
	FELoadCurve(FEModel* fem);
	FELoadCurve(const FELoadCurve& lc);

	void operator = (const FELoadCurve& lc);

	// destructor
	virtual ~FELoadCurve();

	void Serialize(DumpStream& ar);

	bool CopyFrom(FELoadCurve* lc);

	void Add(double time, double value);

	void Clear();

	FEPointFunction& GetFunction() { return m_fnc; }

	void SetInterpolation(FEPointFunction::INTFUNC f);
	void SetExtendMode(FEPointFunction::EXTMODE f);

protected:
	double GetValue(double time) override;

private:
	FEPointFunction	m_fnc;	//!< functin to evaluate

	DECLARE_FECORE_CLASS();
};
