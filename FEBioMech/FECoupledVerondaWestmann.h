#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This material describes a coupled Veronda-Westmann formulation
class FECoupledVerondaWestmann : public FEElasticMaterial
{
public:
	//! constructor
	FECoupledVerondaWestmann(FEModel* pfem) : FEElasticMaterial(pfem){}

public:
	//! data initialization
	void Init();

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent at material point
	tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt);
    
protected:
	double	m_c1;	//!< Veronda-Westmann material parameter c1
	double	m_c2;	//!< Veronda-Westmann material parameter c2
	double	m_k;	//!< bulk-modulus

	DECLARE_PARAMETER_LIST();
};
