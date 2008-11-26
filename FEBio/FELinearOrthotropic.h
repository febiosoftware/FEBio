#pragma once
#include "FEMaterial.h"

class FELinearOrthotropic :	public FEElasticMaterial
{
public:
	double	E1, E2, E3;		// Young's moduli
	double	v12, v23, v31;	// Poisson's ratio
	double	G12, G23, G31;	// Shear moduli

public:
	FELinearOrthotropic() : FEElasticMaterial(FE_LINEAR_ORTHOTROPIC){}

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual void Tangent(double D[6][6], FEMaterialPoint& pt);

	//! data initialization
	void Init();

	//! return bulk modulus
	//! \todo what is the bulk modulus of an orthotropic material?
	double BulkModulus() { return E1/(3.0*(1.0 - 2.0*v12));}

	// declare as registered
	DECLARE_REGISTERED(FELinearOrthotropic);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
