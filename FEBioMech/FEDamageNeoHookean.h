#pragma once
#include "FEElasticMaterial.h"
#include "FEDamageMaterialPoint.h"

//-----------------------------------------------------------------------------
// This material is a first attempt to include damage in hyper-elastic materials.
// It assumes the simple damage model as defined in Simo, CMAME 60 (1987), 153-173

//-----------------------------------------------------------------------------
class FEDamageNeoHookean : public FEElasticMaterial
{
public:
	FEDamageNeoHookean(FEModel* pfem);

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

	double	m_alpha;	//!< damage parameter alpha
	double	m_beta;		//!< damage parameter beta

protected:
	double	m_lam;
	double	m_mu;

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization and checking
	bool Init() override;

	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() override
	{ 
		return new FEDamageMaterialPoint(new FEElasticMaterialPoint);
	}

	// calculate damage reduction factor
	double Damage(FEMaterialPoint& pt);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
