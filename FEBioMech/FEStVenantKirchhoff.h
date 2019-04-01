#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Linear elatic material for large deformations

//! This material can be used when a body undergoes large rotations
//! but small strains.

class FEStVenantKirchhoff : public FEElasticMaterial
{
public:
	FEStVenantKirchhoff(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
