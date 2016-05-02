#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// This is the base class for second-order continuum elastic material
// TODO: This class does not define anything new so we should probably delete it.
class FEElasticMaterial2O : public FEElasticMaterial
{
public:
	FEElasticMaterial2O(FEModel* fem);
};
