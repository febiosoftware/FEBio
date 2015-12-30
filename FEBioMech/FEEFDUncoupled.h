#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for the uncoupled ellipsoidal fiber distribution
class FEEFDUncoupled : public FEUncoupledMaterial
{
public:
	FEEFDUncoupled(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

	//! deviatoric Cauchy stress
	mat3ds DevStress(FEMaterialPoint& pt);

	//! deviatoric spatial tangent
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
public:
	double	m_beta[3];	// power in power-law relation
	double	m_ksi[3];	// coefficient in power-law relation

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
