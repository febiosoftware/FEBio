//
//  FEOsmoticVirialExpansion.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/30/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class that implements osmotic pressure using a virial expansion.
//
class FEOsmoticVirialExpansion : public FEElasticMaterial
{
public:
    FEOsmoticVirialExpansion(FEModel* pfem) : FEElasticMaterial(pfem) { m_c1 = m_c2 = m_c3 = 0; }
    
    //! Initialization routine
    void Init();
    
    //! Returns the Cauchy stress
    virtual mat3ds Stress(FEMaterialPoint& mp);
    
    //! Returs the spatial tangent
    virtual tens4ds Tangent(FEMaterialPoint& mp);
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
    
public:
    double	m_phiwr;	//!< fluid volume fraction in reference configuration
    double	m_cr;		//!< concentration in reference configuration
    double  m_c1;       //!< first virial coefficient
    double  m_c2;       //!< second virial coefficient
    double  m_c3;       //!< third virial coefficient
};
