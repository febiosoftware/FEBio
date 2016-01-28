#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// Material point data for discrete materials.
class FEDiscreteMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy();

	void Init(bool bflag);
};

//-----------------------------------------------------------------------------
//! material class for discrete elements
class FEDiscreteMaterial : public FEMaterial
{
public:
	FEDiscreteMaterial(FEModel* pfem);
	virtual FEMaterialPoint* CreateMaterialPointData();
};
