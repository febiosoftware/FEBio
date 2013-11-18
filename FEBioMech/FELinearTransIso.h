#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a linear transversely isotropic material.
class FELinearTransIso : public FEElasticMaterial
{
public:
	double	E1, E3;		// Young's moduli
	double	v12, v31;	// Poisson's ratio
	double	G23;		// Shear moduli

public:
	FELinearTransIso(FEModel* pfem) : FEElasticMaterial(pfem) {}

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
