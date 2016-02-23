#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEElasticMaterialPoint2O : public FEMaterialPoint
{
public:
	//! constructor
	FEElasticMaterialPoint2O();

	//! initialization
	void init(bool bflag);

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! serialize material point data
	void Serialize(DumpStream& ar);
};

//-----------------------------------------------------------------------------
// This is the base class for second-order continuum elastic material
class FEElasticMaterial2O : public FEElasticMaterial
{
public:
	FEElasticMaterial2O(FEModel* fem);
};
