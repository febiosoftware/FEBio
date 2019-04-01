#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Incompressible Neo-Hookean material

class FEIncompNeoHookean : public FEUncoupledMaterial
{
public:
	FEIncompNeoHookean(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	double	m_G;	//!< Shear modulus

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
