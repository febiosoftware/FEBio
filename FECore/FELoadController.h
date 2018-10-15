#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
// Class that describes a load controller. A load controller can modify the value
// of model parameters during the analysis. 
class FELoadController : public FECoreBase
{
	DECLARE_SUPER_CLASS(FELOADCONTROLLER_ID);

public:
	FELoadController(FEModel* fem);

	//! evaluate the load controller 
	void Evaluate(double time);

	//! return the last calculated value
	double Value() const { return m_value; }

	//! serialization
	void Serialize(DumpStream& ar) override;

protected:
	// This must be implemented by derived classes
	virtual double GetValue(double time) = 0;

private:
	double	m_value;	//!< last calculated value
};
