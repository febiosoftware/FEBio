// FETendonMaterial.h: interface for the FETendonMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FETENDONMATERIAL_H__A76CA7E3_D784_46DC_9466_F46D6DBB8D3C__INCLUDED_)
#define AFX_FETENDONMATERIAL_H__A76CA7E3_D784_46DC_9466_F46D6DBB8D3C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FETransverselyIsotropic.h"

//-----------------------------------------------------------------------------
//! Tendon Material

//! This material uses the constitutive model developed by Blemker et.al. to model
//! tendons
//! Note that lam1 and m_K are inherited

class FETendonMaterial : public FETransverselyIsotropic
{
public:
	FETendonMaterial () {}

public:
	// transverse constants
	double	m_G1;	//!< along-fiber shear modulus
	double	m_G2;	//!< cross-fiber shear modulus

	// along fiber constants
	double	m_L1;	//!< tendon fiber constant L1
	double	m_L2;	//!< tendon fiber constant L2
	
public:
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt);

	// declare as registered
	DECLARE_REGISTERED(FETendonMaterial);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FETENDONMATERIAL_H__A76CA7E3_D784_46DC_9466_F46D6DBB8D3C__INCLUDED_)
