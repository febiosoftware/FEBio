#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a linear transversely isotropic material.
class FELinearTransIso : public FEElasticMaterial
{
public:
	double	E1, E3;		// Young's moduli
	double	v12, v23;	// Poisson's ratio
	double	G12;		// Shear moduli
    double  lam, lamT, lamL;    // Lamé constants
    double  mu, muT;    // Lamé constants

public:
	FELinearTransIso(FEModel* pfem) : FEElasticMaterial(pfem) {}

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
	//! data initialization
	bool Validate();

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
