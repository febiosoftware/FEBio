/*
 *  FESphericalFiberDistribution.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 12/26/11.
 *  Copyright 2011 Columbia University. All rights reserved.
 *
 */

#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for the spherical fiber distribution
//!

class FESphericalFiberDistribution : public FEElasticMaterial
{
public:
	FESphericalFiberDistribution() {m_unstable = true;}
	
	//! Initialization
	void Init();
	
	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp);
	
	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp);
	
	//! returns the bulkmodulus
	double BulkModulus() {return 0;}
	
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_beta;		// power in power-law relation
	double	m_ksi;		// coefficient in power-law relation
	
	static int		m_nres;	// integration rule
	static double	m_cth[];
	static double	m_sth[];
	static double	m_cph[];
	static double	m_sph[];
	static double	m_w[];
};
