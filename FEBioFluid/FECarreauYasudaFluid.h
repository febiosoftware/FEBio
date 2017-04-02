//
//  FECarreauYasudaFluid.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 4/1/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FECarreauYasudaFluid_hpp
#define FECarreauYasudaFluid_hpp

#include "FEViscousFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the viscous stress in a Carreau-Yasuda power-law fluid

class FECarreauYasudaFluid :	public FEViscousFluid
{
public:
    //! constructor
    FECarreauYasudaFluid(FEModel* pfem);
    
    //! viscous stress
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp);
    
    //! tangent of stress with respect to rate of deformation tensor D
    tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp);
    
    //! dynamic viscosity
    double ShearViscosity(FEMaterialPoint& mp);
    
    //! bulk viscosity
    double BulkViscosity(FEMaterialPoint& mp);
    
public:
    double	m_mu0;		//!< shear viscosity at zero shear rate
    double	m_mui;		//!< shear viscosity at infinite shear rate
    double  m_lam;      //!< time constant
    double  m_n;        //!< power-law index
    double  m_a;        //!< exponent
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

#endif /* FECarreauYasudaFluid_hpp */
