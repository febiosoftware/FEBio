/*
 *  FEEllipsoidalFiberDistribution.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *  Copyright 2010 Columbia University. All rights reserved.
 *
 */

#include "FECore/FEMaterial.h"

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

	static int		m_nres;	// integration rule
	static double	m_cth[];
	static double	m_sth[];
	static double	m_cph[];
	static double	m_sph[];
	static double	m_w[];
};
