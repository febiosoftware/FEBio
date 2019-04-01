#pragma once
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
	double			c1;			//!< Mooney-Rivlin coefficient C1
	double			c2;			//!< Mooney-Rivlin coefficient C2
	FEParamVec3		m_fiber;	//!< fiber orientation

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;

protected:
	FEUncoupledFiberExpLinear	m_fib;
	FEActiveFiberContraction*	m_ac;

	// declare parameter list
	DECLARE_FECORE_CLASS();
};
