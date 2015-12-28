// FETCNonlinearOrthotropic.h: interface for the FETCNonlinearOrthotropic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FETCNONLINEARORTHOTROPIC_H__34FDDCF8_45D0_4B57_A0E4_B29EBF0B8411__INCLUDED_)
#define AFX_FETCNONLINEARORTHOTROPIC_H__34FDDCF8_45D0_4B57_A0E4_B29EBF0B8411__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

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
	virtual mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt);

	//! Strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
	//! data initialization and checking
	bool Init();

public:
	double	m_beta[3];
	double	m_ksi[3];

	double m_c1;	//!< Mooney-Rivlin coefficient c1
	double m_c2;	//!< Mooney-Rivlin coefficient c2	

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FETCNONLINEARORTHOTROPIC_H__34FDDCF8_45D0_4B57_A0E4_B29EBF0B8411__INCLUDED_)
