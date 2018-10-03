#pragma once
#include "FEElasticMaterial.h"

class FEGasserOgdenHolzapfel : public FEElasticMaterial
{
public:
	double	m_c;			// neo-Hookean c coefficient
	double	m_k1,m_k2;		// fiber material constants
	double	m_kappa;		// structure coefficient
	double	m_g;			// fiber angle
    double  m_k;            // bulk modulus
		
public:
	FEGasserOgdenHolzapfel(FEModel* pfem) : FEElasticMaterial(pfem) {}
		
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;
		
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
