//
//  FEMembraneMassActionReversible.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 3/4/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FEMembraneMassActionReversible_hpp
#define FEMembraneMassActionReversible_hpp

//-----------------------------------------------------------------------------
//! Law of mass action for reversible membrane reaction
//! (must use effective concentrations).

#include "FEMultiphasic.h"

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

#endif /* FEMembraneMassActionReversible_hpp */
