// FETCNonlinearOrthotropic.h: interface for the FETCNonlinearOrthotropic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FETCNONLINEARORTHOTROPIC_H__34FDDCF8_45D0_4B57_A0E4_B29EBF0B8411__INCLUDED_)
#define AFX_FETCNONLINEARORTHOTROPIC_H__34FDDCF8_45D0_4B57_A0E4_B29EBF0B8411__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! Tension-compression nonlinear orthrotropic

//! This is Gerard's material model for articular cartilage.
//! TODO: make an orthotropic material base class where we 
//! can derive this material from.

class FETCNonlinearOrthotropic : public FEIncompressibleMaterial
{
public:
	FETCNonlinearOrthotropic() {}

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

public:
	double	m_beta[3];
	double	m_ksi[3];

	double m_c1;	//!< Mooney-Rivlin coefficient c1
	double m_c2;	//!< Mooney-Rivlin coefficient c2	

	// declare as registered
	DECLARE_REGISTERED(FETCNonlinearOrthotropic);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FETCNONLINEARORTHOTROPIC_H__34FDDCF8_45D0_4B57_A0E4_B29EBF0B8411__INCLUDED_)
