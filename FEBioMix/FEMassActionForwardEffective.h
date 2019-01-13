//
//  FEMassActionForwardEffective.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 8/23/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEMassActionForwardEffective_hpp
#define FEMassActionForwardEffective_hpp

//-----------------------------------------------------------------------------
//! Law of mass action for forward chemical reaction, using effective concentrations

#include "FEMultiphasic.h"

class FEBIOMIX_API FEMassActionForwardEffective : public FEChemicalReaction
{
public:
    //! constructor
    FEMassActionForwardEffective(FEModel* pfem) : FEChemicalReaction(pfem) {}
    
    //! molar supply at material point
    double ReactionSupply(FEMaterialPoint& pt);
    
    //! tangent of molar supply with strain (J) at material point
    mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt);
    
    //! tangent of molar supply with effective pressure at material point
    double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt);
    
    //! tangent of molar supply with effective concentration at material point
    double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol);
};


#endif /* FEMassActionForwardEffective_hpp */
