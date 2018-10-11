#pragma once
#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
//! Uncoupled formulation of the fiber-exp-linear material for use with uncoupled
//! solid mixtures.
class FEUncoupledFiberExpLinear : public FEElasticFiberMaterialUC
{
public:
	//! Constructor
	FEUncoupledFiberExpLinear(FEModel* pfem);

	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt, const vec3d& n0) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt, const vec3d& n0) override;

	//! calculate deviatoric strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt, const vec3d& n0) override;

public:
	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers
	double	m_lam1;		//!< fiber stretch for straightened fibers

	DECLARE_FECORE_CLASS();
};
