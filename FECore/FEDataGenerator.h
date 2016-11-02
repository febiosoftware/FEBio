#pragma once
#include "vec3d.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEDomain;

//-----------------------------------------------------------------------------
// Data generators are used to generate values of model parameters. 
class FEDataGenerator : public FECoreBase
{
public:
	FEDataGenerator();
	virtual ~FEDataGenerator();

	// this function gives the data generator a chance to initialize itself
	// and check for any input problems.
	virtual bool Init();

	// This function evaluates the variable at all the integration points
	// of the elements in the element set by calling the value function.
	virtual bool Apply(FEDomain* part, const char* szvar);

	// This function does the 
	virtual double value(const vec3d& x) = 0;
};
