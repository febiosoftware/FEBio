//
//  FEFiberPowLinear.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/2/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power law - linear

class FEFiberPowLinear : public FEElasticMaterial
{
public:
    FEFiberPowLinear(FEModel* pfem) : FEElasticMaterial(pfem) { m_thd = 0; m_phd = 90; }
    
    //! Initialization
    void Init();
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp);
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp);
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp);
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
    double  m_beta;     // power law exponent in toe regio
    double	m_thd;		// theta angle for fiber orientation (local coordinates system)
    double	m_phd;		// phi angle for fiber orientation (local coordinates system)
    vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
};
