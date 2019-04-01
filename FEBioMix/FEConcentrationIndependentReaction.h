#pragma once
#include "FEMultiphasic.h"

//-----------------------------------------------------------------------------
//! Concentration-independent forward chemical reaction.

class FEBIOMIX_API FEConcentrationIndependentReaction : public FEChemicalReaction
{
public:
	//! constructor
	FEConcentrationIndependentReaction(FEModel* pfem) : FEChemicalReaction(pfem) {}
		
	//! molar supply at material point
	double ReactionSupply(FEMaterialPoint& pt);
	
	//! tangent of molar supply with strain (J) at material point
	mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective pressure at material point
	double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt);
	
	//! tangent of molar supply with effective concentration at material point
	double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol);
};
