#pragma once
#include "FECoreBase.h"

class FEMaterialPoint;

//---------------------------------------------------------------------------------------
// Base class for evaluating model parameters
class FEValuator : public FECoreBase
{
public:
	FEValuator(FEModel* fem, SUPER_CLASS_ID sid) : FECoreBase(fem, sid) {}
	virtual ~FEValuator() {}
};
