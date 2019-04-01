#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Tension-compression nonlinear orthrotropic

//! This is Gerard's material model for articular cartilage.
//! \todo Make an orthotropic material base class where we 
//!       can derive this material from.

class FETCNonlinearOrthotropic : public FEUncoupledMaterial
{
public:
	FETCNonlinearOrthotropic(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
public:
	double	m_beta[3];
	double	m_ksi[3];

	double m_c1;	//!< Mooney-Rivlin coefficient c1
	double m_c2;	//!< Mooney-Rivlin coefficient c2	

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
