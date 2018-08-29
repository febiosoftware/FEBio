#pragma once

#include "FEUncoupledMaterial.h"
#include "FEDamageMaterialPoint.h"

class FEDamageMooneyRivlin : public FEUncoupledMaterial
{
public:
	FEDamageMooneyRivlin(FEModel* pfem);

public:
	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2

	double	m_beta;		//!< damage parameter beta
	double	m_smin;		//!< damage parameter psi-min
	double	m_smax;		//!< damage parameter psi-max

public:
	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() override { return new FEDamageMaterialPoint(new FEElasticMaterialPoint); }

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
	//! data initialization
	bool Validate() override;

	// calculate damage reduction factor
	double Damage(FEMaterialPoint& pt);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
