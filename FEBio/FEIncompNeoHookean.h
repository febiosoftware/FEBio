// FEIncompNeoHookean.h: interface for the FEIncompNeoHookean class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEINCOMPNEOHOOKEAN_H__8ECAD0BE_54FA_4924_9952_68EA377A8D8E__INCLUDED_)
#define AFX_FEINCOMPNEOHOOKEAN_H__8ECAD0BE_54FA_4924_9952_68EA377A8D8E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! Incompressible Neo-Hookean material

class FEIncompNeoHookean : public FEIncompressibleMaterial
{
public:
	FEIncompNeoHookean() {}

public:
	double	m_G;	//!< Shear modulus

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	// declare as registered
	DECLARE_REGISTERED(FEIncompNeoHookean);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FEINCOMPNEOHOOKEAN_H__8ECAD0BE_54FA_4924_9952_68EA377A8D8E__INCLUDED_)
