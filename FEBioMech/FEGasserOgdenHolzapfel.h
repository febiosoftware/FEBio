#pragma once
#include "FEUncoupledMaterial.h"

class FEGasserOgdenHolzapfel : public FEUncoupledMaterial
{
public:
	double	m_c;			// neo-Hookean c coefficient
	double	m_k1,m_k2;		// fiber material constants
	double	m_kappa;		// structure coefficient
	double	m_g;			// fiber angle
		
public:
	FEGasserOgdenHolzapfel(FEModel* pfem) : FEUncoupledMaterial(pfem) {}
		
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);
	
	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);
		
	//! calculate deviatoric strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt);
    
	//! data initialization
	void Init();
		
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
