#pragma once
#include "FEElasticFiberMaterial.h"

//-----------------------------------------------------------------------------
//! This class represents a fiber material with an exponential toe-region
//! and a linear region.
class FEFiberExpLinear : public FEElasticFiberMaterial
{
public:
	//! constructor
	FEFiberExpLinear(FEModel* pfem);
	
	//! Calculate the fiber stress
	mat3ds Stress(FEMaterialPoint& mp, const vec3d& a0) override;

	//! Calculate the fiber tangent
	tens4ds Tangent(FEMaterialPoint& mp, const vec3d& a0) override;

	//! Calculate the fiber strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) override;

public:
	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers
	double	m_lam1;		//!< fiber stretch for straightened fibers

	DECLARE_FECORE_CLASS();
};
