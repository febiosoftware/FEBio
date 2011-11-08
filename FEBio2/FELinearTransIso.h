#pragma once
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a linear transversely isotropic material.
class FELinearTransIso : public FEElasticMaterial
{
public:
	double	E1, E3;		// Young's moduli
	double	v12, v31;	// Poisson's ratio
	double	G23;		// Shear moduli

public:
	FELinearTransIso() {}

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	//! return bulk modulus
	//! \todo what is the bulk modulus of an trans-iso material?
	double BulkModulus() { return 0;}

	// declare as registered
	DECLARE_REGISTERED(FELinearTransIso);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
