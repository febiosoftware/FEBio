#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// Material point data for discrete materials.
class FECORE_API FEDiscreteMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy();
};

//-----------------------------------------------------------------------------
//! material class for discrete elements
class FECORE_API FEDiscreteMaterial : public FEMaterial
{
public:
	FEDiscreteMaterial(FEModel* pfem);
	virtual FEMaterialPoint* CreateMaterialPointData();
};
