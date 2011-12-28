#pragma once
#include "FEUncoupledMaterial.h"

class FEFungOrthotropic : public FEUncoupledMaterial
{
public:
	double	E1, E2, E3;		// Young's moduli
	double	v12, v23, v31;	// Poisson's ratio
	double	G12, G23, G31;	// Shear moduli
	double	lam[3][3];		// first Lame coefficients
	double	mu[3];			// second Lame coefficients
	double	m_c;			// c coefficient
		
public:
	FEFungOrthotropic() {m_K = 0;}
		
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);
	
	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);
		
	//! data initialization
	void Init();

	//! convert material properties
	void ConvertProperties();

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
