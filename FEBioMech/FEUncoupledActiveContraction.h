#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of an uncoupled solid matrix material.
class FEUncoupledActiveContraction : public FEUncoupledMaterial
{
public:
	//! constructor
	FEUncoupledActiveContraction();

	//! Initialization
	void Init();

	//! deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt);

	//! deviatoric tangent
	tens4ds DevTangent(FEMaterialPoint& pt);

public:
	double	m_Tmax;
	double	m_ca0;
	double	m_camax;
	double	m_beta;
	double	m_l0;
	double	m_refl;

	DECLARE_REGISTERED(FEUncoupledActiveContraction);

	DECLARE_PARAMETER_LIST();
};
