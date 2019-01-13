//
//  FEMembraneReactionRateVoltageGated.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 4/30/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FEMembraneReactionRateVoltageGated_hpp
#define FEMembraneReactionRateVoltageGated_hpp

#include "FEMultiphasic.h"

class FECORE_API FEMembraneReactionRateVoltageGated : public FEMembraneReactionRate
{
public:
    //! constructor
    FEMembraneReactionRateVoltageGated(FEModel* pfem) : FEMembraneReactionRate(pfem) { m_a = m_b = m_c = m_d = 0; m_sol = -1; m_z = 0; }
    
    // initialization
    bool Init() override;
    
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
    double Tangent_ReactionRate_Ce(FEMaterialPoint& pt, const int isol) override;
    double Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol) override;
    
public:
    int     m_sol;      //!< solute id (1-based)
    int     m_z;        //!< charge number of channel ion
    double  m_a;        //!< coefficient
    double  m_b;        //!< coefficient
    double  m_c;        //!< coefficient
    double  m_d;        //!< coefficient
    
    DECLARE_FECORE_CLASS();
};

#endif /* FEMembraneReactionRateVoltageGated_hpp */
