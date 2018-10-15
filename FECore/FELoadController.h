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
	virtual void Evaluate(double time) = 0;

	//! return the last calculated value
	virtual double Value() const = 0;
};
