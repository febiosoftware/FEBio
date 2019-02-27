#pragma once
#include "FECoreBase.h"

class FEMaterialPoint;

class FEModelParam;

//---------------------------------------------------------------------------------------
// Base class for evaluating model parameters
class FECORE_API FEValuator : public FECoreBase
{
public:
	FEValuator(FEModel* fem) : FECoreBase(fem), m_param(nullptr) {}
	virtual ~FEValuator() {}

	void SetModelParam(FEModelParam* p) { m_param = p; }
	FEModelParam* GetModelParam() { return m_param; }

private:
	FEModelParam*	m_param;	//!< the model param that is using this valuator
};
