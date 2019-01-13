//
//  FEMassActionReversibleEffective.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 8/23/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEMassActionReversibleEffective_hpp
#define FEMassActionReversibleEffective_hpp

//-----------------------------------------------------------------------------
//! Law of mass action for reversible chemical reaction
//! using effective concentrations.

#include "FEMultiphasic.h"

class FECORE_API FEMassActionReversibleEffective : public FEChemicalReaction
{
public:
    //! constructor
    FEMassActionReversibleEffective(FEModel* pfem) : FEChemicalReaction(pfem) {}
    
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


#endif /* FEMassActionReversibleEffective_hpp */
