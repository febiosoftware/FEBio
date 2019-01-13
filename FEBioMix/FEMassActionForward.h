/*
 *  FEMassActionForward.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/6/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

//-----------------------------------------------------------------------------
//! Law of mass action for forward chemical reaction.

#include "FEMultiphasic.h"

class FECORE_API FEMassActionForward : public FEChemicalReaction
{
public:
	//! constructor
	FEMassActionForward(FEModel* pfem) : FEChemicalReaction(pfem) {}
		
	//! molar supply at material point
	double ReactionSupply(FEMaterialPoint& pt);
	
	//! tangent of molar supply with strain (J) at material point
	mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective pressure at material point
	double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective concentration at material point
	double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol);
};

