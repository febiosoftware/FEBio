#pragma once
#include "FEUncoupledMaterial.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! Mooney-Rivlin material

class FEMooneyRivlin : public FEUncoupledMaterial
{
public:
	FEMooneyRivlin(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	FEParamDouble	m_c1;	//!< Mooney-Rivlin coefficient C1
	FEParamDouble	m_c2;	//!< Mooney-Rivlin coefficient C2

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
