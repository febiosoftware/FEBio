#pragma once
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a version of orthotropic elasticity that remains
// valid for large deformations
class FEOrthoElastic :	public FEElasticMaterial
{
public:
	double	E1, E2, E3;		// Young's moduli
	double	v12, v23, v31;	// Poisson's ratio
	double	G12, G23, G31;	// Shear moduli
	double	lam[3][3];		// first Lame coefficients
	double	mu[3];			// second Lame coefficients

public:
	FEOrthoElastic() {}

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	//! return bulk modulus
	//! \todo what is the bulk modulus of an orthotropic material?
	double BulkModulus() { return (lam[0][0]+lam[1][1]+lam[2][2]
								   +2*(lam[0][1]+lam[1][2]+lam[0][2]
									   +mu[0]+mu[1]+mu[2]))/9.0;}

	// declare as registered
	DECLARE_REGISTERED(FELinearOrthotropic);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
