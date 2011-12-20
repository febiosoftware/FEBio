#pragma once

#include "FEUncoupledMaterial.h"
#include "FEDamageNeoHookean.h"

class FEDamageMooneyRivlin : public FEUncoupledMaterial
{
public:
	FEDamageMooneyRivlin(void);

public:
	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2

	double	m_beta;		//!< damage parameter beta
	double	m_smin;		//!< damage parameter psi-min
	double	m_smax;		//!< damage parameter psi-max

public:
	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return new FEDamageMaterialPoint(new FEElasticMaterialPoint); }

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	// calculate damage reduction factor
	double Damage(FEMaterialPoint& pt);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
