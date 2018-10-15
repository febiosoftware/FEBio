#pragma once
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This class can be used to define global model data and will be placed in the
//! global date section of the FEModel class
class FECORE_API FEGlobalData : public FECoreBase
{
	DECLARE_SUPER_CLASS(FEGLOBALDATA_ID);

public:
	//! constructor
	FEGlobalData(FEModel* fem);

	// initialization
	virtual bool Init();
};
