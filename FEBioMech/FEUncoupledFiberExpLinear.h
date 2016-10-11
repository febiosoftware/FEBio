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
	mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! calculate deviatoric strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt);

public:
	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers
	double	m_lam1;		//!< fiber stretch for straightened fibers

	DECLARE_PARAMETER_LIST();
};
