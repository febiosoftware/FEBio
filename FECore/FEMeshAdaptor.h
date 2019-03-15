#pragma once
#include "FECoreBase.h"

class FEModel;

class FECORE_API FEMeshAdaptor : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FEMeshAdaptor(FEModel* fem);

	// The mesh adaptor should return true if the mesh remained unchanged
	// otherwise, it should return false
	virtual bool Apply() = 0;
};
