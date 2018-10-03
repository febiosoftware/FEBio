// FEMooneyRivlin.h: interface for the FEMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMOONEYRIVLIN_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_)
#define AFX_FEMOONEYRIVLIN_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEUncoupledMaterial.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! Mooney-Rivlin material

class FEMooneyRivlin : public FEUncoupledMaterial
{
public:
	FEMooneyRivlin(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	FEParamDouble	m_c1;	//!< Mooney-Rivlin coefficient C1
	FEParamDouble	m_c2;	//!< Mooney-Rivlin coefficient C2

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

#endif // !defined(AFX_FEMOONEYRIVLIN_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_)
