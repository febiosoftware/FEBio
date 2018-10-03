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
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
protected:
	double	m_c1;	//!< Veronda-Westmann material parameter c1
	double	m_c2;	//!< Veronda-Westmann material parameter c2
	double	m_k;	//!< bulk-modulus

	DECLARE_FECORE_CLASS();
};
