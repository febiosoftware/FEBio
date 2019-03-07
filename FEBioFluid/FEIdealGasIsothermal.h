//
//  FEIdealGasIsothermal.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEIdealGasIsothermal_hpp
#define FEIdealGasIsothermal_hpp

#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! Ideal gas under isothermal conditions.

class FEBIOFLUID_API FEIdealGasIsothermal : public FEFluid
{
public:
    FEIdealGasIsothermal(FEModel* pfem);
    
public:
    //! initialization
    bool Init() override;
    
    //! elastic pressure
    double Pressure(const double e) override;
    
    //! tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain(FEMaterialPoint& mp) override;
    
    //! 2nd tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp) override;
    
    //! strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! invert pressure-dilatation relation
    double Dilatation(const double p) override;
    
    //! evaluate temperature
    double Temperature(FEMaterialPoint& mp) override;
    
public:
    double      m_M;        //!< moral mass
    double      m_pr;       //!< ambient pressure
    double      m_Tr;       //!< ambient temperature
    double      m_R;        //!< universal gas constant
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

#endif /* FEIdealGasIsothermal_hpp */
