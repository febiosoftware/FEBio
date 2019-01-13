//
//  FESFDParam.h
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 7/16/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#ifndef __FEBioXCode4__FESFDParam__
#define __FEBioXCode4__FESFDParam__

#include "FEBioMech/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for the spherical fiber distribution with
//! fiber modulus dependent on sbm referential density

class FECORE_API FESFDSBM : public FEElasticMaterial
{
public:
	FESFDSBM(FEModel* pfem) : FEElasticMaterial(pfem) { m_alpha = 0;}
	
	//! Initialization
	bool Init() override;
	
	//! Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;
	
	// Spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	// Strain energy density
	virtual double StrainEnergyDensity(FEMaterialPoint& mp) override;
	
	//! return fiber modulus
	double FiberModulus(double rhor) { return m_ksi0*pow(rhor/m_rho0, m_g);}
	
	// declare the parameter list
	DECLARE_FECORE_CLASS();
	
public:
	double	m_alpha;	// coefficient of exponential argument
	double	m_beta;		// power in power-law relation
	double	m_ksi0;		// ksi = ksi0*(rhor/rho0)^gamma
    double  m_rho0;     // rho0
    double  m_g;        // gamma
	int		m_sbm;      //!< global id of solid-bound molecule
	int		m_lsbm;     //!< local id of solid-bound molecule
	
	static int		m_nres;	// integration rule
	static double	m_cth[];
	static double	m_sth[];
	static double	m_cph[];
	static double	m_sph[];
	static double	m_w[];
};

#endif /* defined(__FEBioXCode4__FESFDParam__) */
