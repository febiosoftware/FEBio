// FEMooneyRivlin.h: interface for the FEMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Coupled transversely-isotropic Mooney-Rivlin material
//!
class FECoupledTransIsoMooneyRivlin : public FEElasticMaterial
{
public:
	FECoupledTransIsoMooneyRivlin(FEModel* pfem) : FEElasticMaterial(pfem) {}

public:
	double	m_c1;	//!< Mooney-Rivlin coefficient C1
	double	m_c2;	//!< Mooney-Rivlin coefficient C2
	double	m_c3;	//!< fiber stress scale factor
	double	m_c4;	//!< exponential scale factor
	double	m_c5;	//!< slope of linear stress region
	double	m_flam;	//!< fiber stretch at which fibers are straight
	double	m_K;	//!< "bulk"-modulus

private:
	//! calculate deviatoric stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
public:
	tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
