// FEIncompNeoHookean.h: interface for the FEIncompNeoHookean class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEINCOMPNEOHOOKEAN_H__8ECAD0BE_54FA_4924_9952_68EA377A8D8E__INCLUDED_)
#define AFX_FEINCOMPNEOHOOKEAN_H__8ECAD0BE_54FA_4924_9952_68EA377A8D8E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Incompressible Neo-Hookean material

class FEIncompNeoHookean : public FEUncoupledMaterial
{
public:
	FEIncompNeoHookean() {}

public:
	double	m_G;	//!< Shear modulus

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! initialization
	void Init();

	// declare as registered
	DECLARE_REGISTERED(FEIncompNeoHookean);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FEINCOMPNEOHOOKEAN_H__8ECAD0BE_54FA_4924_9952_68EA377A8D8E__INCLUDED_)
