#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! This material was added by Shawn Reese

class FENeoHookeanTransIso : public FEElasticMaterial
{
public:
	FENeoHookeanTransIso(void) : FEElasticMaterial(FE_NEOHOOKEAN_TRANSISO){}

public:
	double	m_Ep;	//!< Young's modulus
	double	m_Ez;	//!< Young's modulus
	double	m_vz;	//!< Poisson's ratio
	double	m_vp;	//!< Poisson's ratio
	double	m_gz;	//!< shear modulus

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual void Tangent(double D[6][6], FEMaterialPoint& pt);

	double BulkModulus() { return 0; }

	// declare as registered
	DECLARE_REGISTERED(FENeoHookeanTransIso);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
