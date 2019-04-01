#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Arruda-Boyce material

class FEArrudaBoyce : public FEUncoupledMaterial
{
public:
	FEArrudaBoyce(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	FEParamDouble	m_mu;	//!< shear modulus
	double	m_N;	//!< Nr of links in chain

public:
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density
	virtual double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
