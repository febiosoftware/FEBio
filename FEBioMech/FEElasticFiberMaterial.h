#pragma once
#include "FEElasticMaterial.h"
#include <FECore/FEVectorGenerator.h>

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterial : public FEMaterial
{
public:
    FEElasticFiberMaterial(FEModel* pfem);

	// calculate stress in fiber direction a0
	virtual mat3ds Stress(FEMaterialPoint& mp, const vec3d& a0) = 0;

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp, const vec3d& a0) = 0;

	//! Strain energy density
	virtual double StrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) = 0;
};
