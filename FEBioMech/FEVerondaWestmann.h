#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//!  Veronda-Westmann material model

class FEVerondaWestmann : public FEUncoupledMaterial
{
public:
	FEVerondaWestmann(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	FEParamDouble	m_c1;	//!< Veronda-Westmann coefficient C1;
	FEParamDouble	m_c2;	//!< Veronda-Westmann coefficient C2;

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
