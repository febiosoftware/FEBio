#pragma once
#include "FEMultiphasic.h"

//-----------------------------------------------------------------------------
//! Law of mass action for forward chemical reaction, using effective concentrations
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
