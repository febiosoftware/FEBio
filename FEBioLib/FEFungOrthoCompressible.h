#pragma once
#include "FECore/FEElasticMaterial.h"

class FEFungOrthoCompressible : public FEElasticMaterial
{
public:
	double	E1, E2, E3;		// Young's moduli
	double	v12, v23, v31;	// Poisson's ratio
	double	G12, G23, G31;	// Shear moduli
	double	lam[3][3];		// first Lame coefficients
	double	mu[3];			// second Lame coefficients
	double	m_c;			// c coefficient
	double	m_k;			// bulk modulus
	
public:
	FEFungOrthoCompressible() {}
	
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);

	//! convert material properties
	void ConvertProperties();
	
	//! data initialization
	void Init();
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
