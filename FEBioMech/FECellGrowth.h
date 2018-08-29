/*
 *  FECellGrowth.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 7/8/11.
 *  Copyright 2011 Columbia University. All rights reserved.
 *
 */

#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class that implements the equilibrium of a perfect osmometer.
//
class FECellGrowth : public FEElasticMaterial
{
public:
	//! When used on its own (not in a solid mixture), this materials
	//! is intrinsically unstable
	FECellGrowth(FEModel* pfem) : FEElasticMaterial(pfem) { m_Rgas = 0; m_Tabs = 0; }
	
	//! Initialization routine
	bool Init() override;
	
	//! Returns the Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;
	
	//! Returs the spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;
	
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_phir;		//!< intracellular solid volume fraction normalized to reference configuration
	double	m_cr;		//!< intracellular osmolarity normalized to reference configuration
	double	m_ce;		//!< extracellular osmolarity
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
};
