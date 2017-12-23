// FETransIsoMooneyRivlin.h: interface for the FETransIsoMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FETRANSISOMOONEYRIVLIN_H__E918D89B_4CCD_44B9_9731_19CEC4EDF406__INCLUDED_)
#define AFX_FETRANSISOMOONEYRIVLIN_H__E918D89B_4CCD_44B9_9731_19CEC4EDF406__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEUncoupledMaterial.h"
#include "FEUncoupledFiberExpLinear.h"
#include "FEActiveFiberContraction.h"

//-----------------------------------------------------------------------------
//! Transversely Isotropic Mooney-Rivlin material

//! This material has an isotopric Mooney-Rivlin basis and single preferred
//! fiber direction.

class FETransIsoMooneyRivlin: public FEUncoupledMaterial
{
public:
	FETransIsoMooneyRivlin(FEModel* pfem);

public:
	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;

	//! Create material point data
	FEMaterialPoint* CreateMaterialPointData() override;
    
protected:
	FEUncoupledFiberExpLinear				m_fib;
	FEPropertyT<FEActiveFiberContraction>	m_ac;

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FETRANSISOMOONEYRIVLIN_H__E918D89B_4CCD_44B9_9731_19CEC4EDF406__INCLUDED_)
