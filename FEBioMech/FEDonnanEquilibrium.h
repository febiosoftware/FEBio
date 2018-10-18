/*
 *  FEDonnanEquilibrium.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *
 */
#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class that implements Donnan equilibrium. 
//! When used on its own (not in a solid mixture), this materials
//! is intrinsically unstable
class FEDonnanEquilibrium : public FEElasticMaterial
{
public:
	// constructor
	FEDonnanEquilibrium(FEModel* pfem);
	
	//! Initialization routine
	bool Init() override;

	//! Returns the Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;

	//! Returs the spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
	
public:
	double	m_phiwr;	//!< fluid volume fraction in reference configuration
    double	m_phisr;	//!< referential solid volume fraction (may evolve with time)
	double	m_cFr;		//!< fixed charge density in reference configuration
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
	double	m_bosm;		//!< bath osmolarity
    double  m_Phi;      //!< osmotic coefficient
    bool    m_bnew;     //!< flag for using old or new method
    bool    m_binit;    //!< initialization flag
};
