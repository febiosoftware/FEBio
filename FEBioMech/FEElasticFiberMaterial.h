#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterial : public FEElasticMaterial
{
public:
    FEElasticFiberMaterial(FEModel* pfem);

	// calculate stress in fiber direction a0
	virtual mat3ds Stress(FEMaterialPoint& mp, const vec3d& a0) = 0;

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp, const vec3d& a0) = 0;

	//! Strain energy density
	virtual double StrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) = 0;

private:
	// These are made private since fiber materials should implement the functions above instead. 
	// The functions can still be reached when a fiber material is used in an elastic mixture. 
	// In those cases the fiber vector is taken from the first column of Q. 
	mat3ds Stress(FEMaterialPoint& mp) final;
	tens4ds Tangent(FEMaterialPoint& mp) final;
	double StrainEnergyDensity(FEMaterialPoint& mp) final;
};
