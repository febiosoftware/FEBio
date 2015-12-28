//
//  FEEFD.h
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 5/18/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#ifndef __FEBioXCode4__FEEFD__
#define __FEBioXCode4__FEEFD__

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for the ellipsoidal fiber distribution
//!

class FEEFD : public FEElasticMaterial
{
public:
	FEEFD(FEModel* pfem) : FEElasticMaterial(pfem) {}
	
	//! Initialization
	bool Init();
    
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);
	
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_beta[3];	// power in power-law relation
	double	m_ksi[3];	// coefficient in power-law relation
};


#endif /* defined(__FEBioXCode4__FEEFD__) */
