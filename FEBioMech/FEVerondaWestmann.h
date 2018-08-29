// FEVerondaWestmann.h: interface for the FEVerondaWestmann class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEVERONDAWESTMANN_H__0BA871E2_75AF_426D_BB95_B09FECFB5C9A__INCLUDED_)
#define AFX_FEVERONDAWESTMANN_H__0BA871E2_75AF_426D_BB95_B09FECFB5C9A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//!  Veronda-Westmann material model

class FEVerondaWestmann : public FEUncoupledMaterial
{
public:
	FEVerondaWestmann(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	double	m_c1;	//!< Veronda-Westmann coefficient C1;
	double	m_c2;	//!< Veronda-Westmann coefficient C2;

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FEVERONDAWESTMANN_H__0BA871E2_75AF_426D_BB95_B09FECFB5C9A__INCLUDED_)
