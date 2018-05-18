#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This class can be used to define global model data and will be placed in the
//! global date section of the FEModel class
class FECORE_API FEGlobalData : public FECoreBase
{
public:
	//! constructor
	FEGlobalData(FEModel* pfem);

	// initialization
	virtual bool Init();

	// get the FEModel
	FEModel* GetFEModel() { return m_pfem; }

private:
	FEModel*	m_pfem;
};
