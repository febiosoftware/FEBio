// FEVerondaWestmann.h: interface for the FEVerondaWestmann class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEVERONDAWESTMANN_H__0BA871E2_75AF_426D_BB95_B09FECFB5C9A__INCLUDED_)
#define AFX_FEVERONDAWESTMANN_H__0BA871E2_75AF_426D_BB95_B09FECFB5C9A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//!  Veronda-Westmann material model

class FEVerondaWestmann : public FEIncompressibleMaterial
{
public:
	FEVerondaWestmann() {}

public:
	double	m_c1;	//!< Veronda-Westmann coefficient C1;
	double	m_c2;	//!< Veronda-Westmann coefficient C2;

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! initialize
	void Init();

	// declare as registered
	DECLARE_REGISTERED(FEVerondaWestmann);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FEVERONDAWESTMANN_H__0BA871E2_75AF_426D_BB95_B09FECFB5C9A__INCLUDED_)
