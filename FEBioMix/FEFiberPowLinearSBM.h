//
//  FEFiberPowLinearSBM.h
//  FEBioMix
//
//  Created by Gerard Ateshian on 5/24/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMix__FEFiberPowLinearSBM__
#define __FEBioMix__FEFiberPowLinearSBM__

#include "FEBioMech/FEElasticMaterial.h"
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power law - linear
//! Fiber modulus depends on SBM content

class FEBIOMIX_API FEFiberPowLinearSBM : public FEElasticMaterial
{
public:
    FEFiberPowLinearSBM(FEModel* pfem) : FEElasticMaterial(pfem) { m_thd = 0; m_phd = 90; m_sbm = 0; }
    
    //! Initialization
    bool Init() override;
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp) override;
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp) override;
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! return fiber modulus
    double FiberModulus(double rhor) { return m_E0*pow(rhor/m_rho0, m_g);}
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
    
public:
    double	m_E0;		// fiber modulus E = E0*(rhor/rho0)^gamma
    double  m_lam0;     // stretch ratio at end of toe region
    double  m_beta;     // power law exponent in toe region
    double  m_rho0;     // rho0
    double  m_g;        // gamma
    int		m_sbm;      //!< global id of solid-bound molecule
    int		m_lsbm;     //!< local id of solid-bound molecule
    double	m_thd;		// theta angle for fiber orientation (local coordinates system)
    double	m_phd;		// phi angle for fiber orientation (local coordinates system)
    vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
};

#endif /* defined(__FEBioMix__FEFiberPowLinearSBM__) */
