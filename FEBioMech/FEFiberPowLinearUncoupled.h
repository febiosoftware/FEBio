//
//  FEFiberPowLinearUncoupled.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 5/6/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power-law linear response (uncoupled)

class FEFiberPowLinearUncoupled : public FEUncoupledMaterial
{
public:
    FEFiberPowLinearUncoupled(FEModel* pfem) : FEUncoupledMaterial(pfem) { m_thd = 0; m_phd = 90; }
    
    //! Initialization
    bool Init();
    
    //! Cauchy stress
    virtual mat3ds DevStress(FEMaterialPoint& mp);
    
    // Spatial tangent
    virtual tens4ds DevTangent(FEMaterialPoint& mp);
    
    //! Strain energy density
    virtual double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
    double  m_I0;       // m_lam0^2
    double  m_ksi;      // power law coefficient in toe region
    double  m_beta;     // power law exponent in toe region
    double  m_b;        // coefficient in linear region
    double	m_thd;		// theta angle for fiber orientation (local coordinates system)
    double	m_phd;		// phi angle for fiber orientation (local coordinates system)
    vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
};
