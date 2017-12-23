#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of an uncoupled solid matrix material.
class FEUncoupledActiveContraction : public FEUncoupledMaterial
{
public:
	//! constructor
	FEUncoupledActiveContraction(FEModel* pfem);

	//! deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! deviatoric tangent
	tens4ds DevTangent(FEMaterialPoint& pt) override;

public:
	double	m_Tmax;
	double	m_ca0;
	double	m_camax;
	double	m_beta;
	double	m_l0;
	double	m_refl;

	DECLARE_PARAMETER_LIST();
};
