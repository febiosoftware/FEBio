//
//  FEFiberExpPowSBM.h
//  FEBioMix
//
//  Created by Gerard Ateshian on 5/24/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMix__FEFiberExpPowSBM__
#define __FEBioMix__FEFiberExpPowSBM__

#include "FEBioMech/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Exponential-power law
//! Fiber modulus depends on SBM content

class FEFiberExpPowSBM : public FEElasticMaterial
{
public:
    FEFiberExpPowSBM(FEModel* pfem) : FEElasticMaterial(pfem) { m_thd = 0; m_phd = 90; m_sbm = 0; }
    
    //! Initialization
    bool Init() override;
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp) override;
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp) override;
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! return fiber modulus
    double FiberModulus(double rhor) { return m_ksi0*pow(rhor/m_rho0, m_g);}
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
    
public:
    double	m_alpha;	// coefficient of (In-1) in exponential
    double	m_beta;		// power of (In-1) in exponential
    double	m_ksi0;		// fiber modulus ksi = ksi0*(rhor/rho0)^gamma
    double  m_rho0;     // rho0
    double  m_g;        // gamma
    int		m_sbm;      //!< global id of solid-bound molecule
    int		m_lsbm;     //!< local id of solid-bound molecule
    double	m_thd;		// theta angle for fiber orientation (local coordinates system)
    double	m_phd;		// phi angle for fiber orientation (local coordinates system)
    vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
};

#endif /* defined(__FEBioMix__FEFiberExpPowSBM__) */
