/*
 *  FEDonnanEquilibrium.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *
 */
#pragma once
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class that implements Donnan equilibrium. 
//
class FEDonnanEquilibrium : public FEElasticMaterial
{
public:
	//! When used on its own (not in a solid mixture), this materials
	//! is intrinsically unstable
	FEDonnanEquilibrium() {m_unstable = true; m_Rgas = 0; m_Tabs = 0; }
	
	//! Initialization routine
	void Init();

	//! Returns the Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp);

	//! Returs the spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp);

	//! returns the bulkmodulus
	double BulkModulus();

	// declare as registered
	DECLARE_REGISTERED(FEDonnanEquilibrium);
	
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_phiwr;	//!< fluid volume fraction in reference configuration
	double	m_cFr;		//!< fixed charge density in reference configuration
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
	double	m_bosm;		//!< bath osmolarity
};
