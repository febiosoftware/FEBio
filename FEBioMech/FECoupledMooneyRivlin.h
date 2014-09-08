#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This is a coupled formulation for the Mooney-Rivlin material.
class FECoupledMooneyRivlin : public FEElasticMaterial
{
public:
	FECoupledMooneyRivlin(FEModel* pfem) : FEElasticMaterial(pfem){}

protected:
	double	m_c1;	//!< Mooney-Rivlin parameter c1
	double	m_c2;	//!< Mooney-Rivlin parameter c2
	double	m_K;	//!< bulk modulus

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent at material point
	tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt);
    
	//! data initialization
	void Init();

	DECLARE_PARAMETER_LIST();
};
