// FEMooneyRivlin.h: interface for the FEMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEARRUDABOYCE_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_)
#define AFX_FEARRUDABOYCE_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Arruda-Boyce material

class FEArrudaBoyce : public FEUncoupledMaterial
{
public:
	FEArrudaBoyce(FEModel* pfem) : FEUncoupledMaterial(pfem) {}

public:
	FEParamDouble	m_mu;	//!< shear modulus
	double	m_N;	//!< Nr of links in chain

public:
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density
	virtual double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

#endif // !defined(AFX_FEARRUDABOYCE_H__D75E80FD_0E25_4A12_9539_044C9DC4CB41__INCLUDED_)
