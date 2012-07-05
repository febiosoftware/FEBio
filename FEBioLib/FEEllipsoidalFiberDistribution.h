#pragma once
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for the ellipsoidal fiber distribution
//!

class FEEllipsoidalFiberDistribution : public FEElasticMaterial
{
public:
	FEEllipsoidalFiberDistribution() {m_unstable = true;}
	
	//! Initialization
	void Init();

	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp);

	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp);
	
	//! returns the bulkmodulus
	double BulkModulus() {return 0;}
	
	// declare as registered
	DECLARE_REGISTERED(FEEllipsoidalFiberDistribution);
	
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_beta[3];	// power in power-law relation
	double	m_ksi[3];	// coefficient in power-law relation
};
