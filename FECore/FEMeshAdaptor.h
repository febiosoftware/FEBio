#pragma once
#include "FECoreBase.h"

class FEModel;

class FECORE_API FEMeshAdaptor : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FEMeshAdaptor(FEModel* fem);

	// The mesh adaptor should return true if the mesh remained unchanged
	// otherwise, it should return false.
	// iteration is the iteration number of the mesh adaptation loop
	virtual bool Apply(int iteration) = 0;
};
