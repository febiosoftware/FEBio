//
//  FEMembraneReactionRateIonChannel.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 4/6/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FEMembraneReactionRateIonChannel_hpp
#define FEMembraneReactionRateIonChannel_hpp

#include "FEMultiphasic.h"

class FEMembraneReactionRateIonChannel : public FEMembraneReactionRate
{
public:
    //! constructor
    FEMembraneReactionRateIonChannel(FEModel* pfem) : FEMembraneReactionRate(pfem) { m_g = 0; m_sol = -1; }
    
    //! reaction rate at material point
    double ReactionRate(FEMaterialPoint& pt) override;
    
    //! tangent of reaction rate with strain at material point
    double Tangent_ReactionRate_Strain(FEMaterialPoint& pt) override { return 0; }
    
    //! tangent of reaction rate with effective fluid pressure at material point
    double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) override {return 0; }
    double Tangent_ReactionRate_Pe(FEMaterialPoint& pt) override { return 0; }
    double Tangent_ReactionRate_Pi(FEMaterialPoint& pt) override { return 0; }
    
    //! tangent of reaction rate with effective solute concentration at material point
    double Tangent_ReactionRate_Concentration(FEMaterialPoint& pt, const int isol) override {return 0; }
    double Tangent_ReactionRate_Ce(FEMaterialPoint& pt, const int isol) override { return 0; }
    double Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol) override;
    
public:
    int     m_sol;      //!< solute id (1-based)
    double  m_g;        //!< channel conductance
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEMembraneReactionRateIonChannel_hpp */
