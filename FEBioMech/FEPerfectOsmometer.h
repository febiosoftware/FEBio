/*
 *  FEPerfectOsmometer.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 7/14/10.
 *
 */
#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class that implements the equilibrium of a perfect osmometer.
//
class FEPerfectOsmometer : public FEElasticMaterial
{
public:
	//! When used on its own (not in a solid mixture), this materials
	//! is intrinsically unstable
	FEPerfectOsmometer(FEModel* pfem) : FEElasticMaterial(pfem) {m_unstable = false; m_Rgas = 0; m_Tabs = 0; }
	
	//! Initialization routine
	void Init();

	//! Returns the Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp);

	//! Returs the spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_phiwr;	//!< fluid volume fraction in reference configuration
	double	m_iosm;		//!< internal osmolarity in reference configuration
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
	double	m_bosm;		//!< bath osmolarity
};
