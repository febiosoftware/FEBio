/*
 *  FEReactionRateHuiskes.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 5/15/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "FEMultiphasic.h"

class FEReactionRateHuiskes : public FEReactionRate
{
public:
	//! constructor
	FEReactionRateHuiskes(FEModel* pfem) : FEReactionRate(pfem) { m_B = m_psi0 = 0; }
	
	//! reaction rate at material point
	double ReactionRate(FEMaterialPoint& pt) override;
	
	//! tangent of reaction rate with strain at material point
	mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt) override;
	
	//! tangent of reaction rate with effective fluid pressure at material point
	double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) override;
	
public:
	double	m_B;					//!< mass supply coefficient
	double	m_psi0;					//!< specific strain energy at homeostasis
	
	DECLARE_FECORE_CLASS();	
};
