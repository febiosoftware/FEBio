#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Tendon Material

//! This material uses the constitutive model developed by Blemker et.al. to model tendons
class FETendonMaterial : public FEUncoupledMaterial
{
public:
	FETendonMaterial(FEModel* pfem);

public:
	// transverse constants
	double	m_G1;	//!< along-fiber shear modulus
	double	m_G2;	//!< cross-fiber shear modulus

	// along fiber constants
	double	m_L1;	//!< tendon fiber constant L1
	double	m_L2;	//!< tendon fiber constant L2

	double	m_lam1;	//!< max exponential fiber stretch
	
public:
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
