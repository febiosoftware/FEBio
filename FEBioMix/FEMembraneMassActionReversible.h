#pragma once
#include "FEMultiphasic.h"

//-----------------------------------------------------------------------------
//! Law of mass action for reversible membrane reaction
//! (must use effective concentrations).

class FEBIOMIX_API FEMembraneMassActionReversible : public FEMembraneReaction
{
public:
    //! constructor
    FEMembraneMassActionReversible(FEModel* pfem) : FEMembraneReaction(pfem) {}
    
    //! molar supply at material point
    double ReactionSupply(FEMaterialPoint& pt);
    
    //! tangent of molar supply with strain at material point
    double Tangent_ReactionSupply_Strain(FEMaterialPoint& pt);
    
    //! tangent of molar supply with effective pressure at material point
    double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt);
    double Tangent_ReactionSupply_Pi(FEMaterialPoint& pt);
    double Tangent_ReactionSupply_Pe(FEMaterialPoint& pt);

    //! tangent of molar supply with effective concentration at material point
    double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol);
    double Tangent_ReactionSupply_Ci(FEMaterialPoint& pt, const int sol);
    double Tangent_ReactionSupply_Ce(FEMaterialPoint& pt, const int sol);
    
    //! molar supply at material point
    double FwdReactionSupply(FEMaterialPoint& pt);
    
    //! molar supply at material point
    double RevReactionSupply(FEMaterialPoint& pt);
};
