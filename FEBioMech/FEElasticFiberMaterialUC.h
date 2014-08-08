//
//  FEElasticFiberMaterialUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEElasticFiberMaterialUC__
#define __FEBioMech__FEElasticFiberMaterialUC__

#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterialUC : public FEUncoupledMaterial
{
public:
    FEElasticFiberMaterialUC(FEModel* pfem) : FEUncoupledMaterial(pfem) {}
    
    void SetFiberDirection(FEMaterialPoint& mp, const vec3d n0);
};

//-----------------------------------------------------------------------------
//! Exponential-power law

class FEFiberExponentialPowerUC : public FEElasticFiberMaterialUC
{
public:
	FEFiberExponentialPowerUC(FEModel* pfem) : FEElasticFiberMaterialUC(pfem) {}
	
	//! Initialization
	void Init();
	
	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp);
	
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp);
	
public:
	double	m_alpha;	// coefficient of (In-1) in exponential
	double	m_beta;		// power of (In-1) in exponential
	double	m_ksi;		// fiber modulus
    double  m_mu;       // shear modulus
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Neo-Hookean law

class FEFiberNHUC : public FEElasticFiberMaterialUC
{
public:
	FEFiberNHUC(FEModel* pfem) : FEElasticFiberMaterialUC(pfem) {}
	
	//! Initialization
	void Init();
	
	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp);
	
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp);
	
public:
	double	m_mu;       // shear modulus
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEElasticFiberMaterialUC__) */
