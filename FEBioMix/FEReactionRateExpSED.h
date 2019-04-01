#pragma once
#include "FEMultiphasic.h"

class FEBIOMIX_API FEReactionRateExpSED : public FEReactionRate
{
public:
    //! constructor
    FEReactionRateExpSED(FEModel* pfem) : FEReactionRate(pfem) { m_B = m_Psi0 = 0; }
    
    //! reaction rate at material point
    double ReactionRate(FEMaterialPoint& pt) override;
    
    //! tangent of reaction rate with strain at material point
    mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt) override;
    
    //! tangent of reaction rate with effective fluid pressure at material point
    double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) override;
    
public:
    double	m_B;					//!< mass supply coefficient
    double	m_Psi0;					//!< scaling strain energy density
    
    DECLARE_FECORE_CLASS();	
};
