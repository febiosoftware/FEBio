//
//  FEReactionRateExpSED.h
//  FEBioMix
//
//  Created by Gerard Ateshian on 2/3/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMix__FEReactionRateExpSED__
#define __FEBioMix__FEReactionRateExpSED__

#include "FEMultiphasic.h"

class FEReactionRateExpSED : public FEReactionRate
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

#endif /* defined(__FEBioMix__FEReactionRateExpSED__) */
