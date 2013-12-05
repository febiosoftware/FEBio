//
//  FEConcentrationIndependentReaction.h
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 7/21/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#ifndef __FEBioXCode4__FEConcentrationIndependentReaction__
#define __FEBioXCode4__FEConcentrationIndependentReaction__

//-----------------------------------------------------------------------------
//! Concentration-independent forward chemical reaction.

#include "FEMultiphasic.h"

class FEConcentrationIndependentReaction : public FEChemicalReaction
{
public:
	//! constructor
	FEConcentrationIndependentReaction(FEModel* pfem) : FEChemicalReaction(pfem) {}
	
	//! data initialization and checking
	void Init() {FEChemicalReaction::Init(); }
	
	//! molar supply at material point
	double ReactionSupply(FEMaterialPoint& pt);
	
	//! tangent of molar supply with strain (J) at material point
	mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective pressure at material point
	double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective concentration at material point
	double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol);
};

#endif /* defined(__FEBioXCode4__FEConcentrationIndependentReaction__) */
