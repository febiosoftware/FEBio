//
//  FEIdealGasIsentropic.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/11/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FEIdealGasIsentropic_hpp
#define FEIdealGasIsentropic_hpp

#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! Base class for fluid materials.

class FEIdealGasIsentropic : public FEFluid
{
public:
    FEIdealGasIsentropic(FEModel* pfem);
    
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
    
private: // material properties
    FEViscousFluid* m_pViscous; //!< pointer to viscous part of fluid material
    
public:
    double      m_gamma;    //!< ratio of specific heats (constant pressure/constant volume)
    double      m_M;        //!< moral mass
    double      m_pr;       //!< ambient pressure
    double      m_Tr;       //!< ambient temperature
    double      m_R;        //!< universal gas constant
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

#endif /* FEIdealGasIsentropic_hpp */
