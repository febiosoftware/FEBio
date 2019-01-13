/*
 *  FEMassActionReversible.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/7/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

//-----------------------------------------------------------------------------
//! Law of mass action for reversible chemical reaction.

#include "FEMultiphasic.h"

class FECORE_API FEMassActionReversible : public FEChemicalReaction
{
public:
	//! constructor
	FEMassActionReversible(FEModel* pfem) : FEChemicalReaction(pfem) {}
		
	//! molar supply at material point
	double ReactionSupply(FEMaterialPoint& pt);
	
	//! tangent of molar supply with strain (J) at material point
	mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective pressure at material point
	double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective concentration at material point
	double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol);
	
	//! molar supply at material point
	double FwdReactionSupply(FEMaterialPoint& pt);
	
	//! molar supply at material point
	double RevReactionSupply(FEMaterialPoint& pt);	
};

