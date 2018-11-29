// FEMuscleMaterial.h: interface for the FEMuscleMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMUSCLEMATERIAL_H__528059E0_10E8_49A1_8168_DB5EFBEB2A93__INCLUDED_)
#define AFX_FEMUSCLEMATERIAL_H__528059E0_10E8_49A1_8168_DB5EFBEB2A93__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEUncoupledMaterial.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! Muscle Material

//! This material uses the constitutive model developed by Blemker et.al. to model
//! muscles which undergo active contraction

class FEMuscleMaterial: public FEUncoupledMaterial
{
public:
	FEMuscleMaterial(FEModel* pfem);

public:
	// transverse constants
	double m_G1; //!< along-fiber shear modulus
	double m_G2; //!< cross-fiber shear modulus
	double m_G3; //!< new term

	// along fiber constants
	double	m_P1; //!< muscle fiber constant P1
	double	m_P2; //!< muscle fiber constant P2
		
	double	m_Lofl;  //!< optimal sarcomere length
	double	m_smax;  //!< maximum isometric stretch
	double	m_lam1;
	double	m_alpha;	//!< activation parameter

	FEParamVec3		m_fiber;

public:
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	// declare the material parameters
	DECLARE_FECORE_CLASS();
};


#endif // !defined(AFX_FEMUSCLEMATERIAL_H__528059E0_10E8_49A1_8168_DB5EFBEB2A93__INCLUDED_)
