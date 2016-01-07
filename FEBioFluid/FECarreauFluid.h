//
//  FECarreauFluid.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/29/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#pragma once
#include "FEViscousFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the viscous stress in a Carreau power-law fluid

class FECarreauFluid :	public FEViscousFluid
{
public:
    //! constructor
    FECarreauFluid(FEModel* pfem);
    
    //! viscous stress
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp);
    
    //! tangent of stress with respect to rate of deformation tensor D
    tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp);
    
    //! dynamic viscosity
    double DynamicViscosity(FEMaterialPoint& mp);
    
public:
    double	m_mu0;		//!< shear viscosity at zero shear rate
    double	m_mui;		//!< shear viscosity at infinite shear rate
    double  m_lam;      //!< time constant
    double  m_n;        //!< power-law index
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};
