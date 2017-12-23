#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a linear orthotropic material.
class FELinearOrthotropic :	public FEElasticMaterial
{
public:
	double	E1, E2, E3;		// Young's moduli
	double	v12, v23, v31;	// Poisson's ratio
	double	G12, G23, G31;	// Shear moduli
	double	lam[3][3];		// first Lame coefficients
	double	mu[3];			// second Lame coefficients

public:
	FELinearOrthotropic(FEModel* pfem) : FEElasticMaterial(pfem) {}

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization
	bool Validate() override;

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
