// FEMooneyRivlin.h: interface for the FEMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Coupled transversely-isotropic Mooney-Rivlin material
//!
class FECoupledTransIsoMooneyRivlin : public FEElasticMaterial
{
public:
	FECoupledTransIsoMooneyRivlin(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	double	m_c1;	//!< Mooney-Rivlin coefficient C1
	double	m_c2;	//!< Mooney-Rivlin coefficient C2
	double	m_c3;	//!< fiber stress scale factor
	double	m_c4;	//!< exponential scale factor
	double	m_c5;	//!< slope of linear stress region
	double	m_flam;	//!< fiber stretch at which fibers are straight
	double	m_K;	//!< "bulk"-modulus

public:
	//! calculate deviatoric stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
