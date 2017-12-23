//
//  FECubicCLE.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/26/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a cubic conewise linear elastic (CLE) material (Curnier et al. 1995 J Elasticity).
class FECubicCLE :	public FEElasticMaterial
{
public:
    double	m_lp1;      // diagonal first lamé constants (tension)
    double	m_lm1;      // diagonal first lamé constants (compression)
    double	m_l2;       // off-diagonal first lamé constants
    double	m_mu;       // shear moduli
    
public:
    FECubicCLE(FEModel* pfem) : FEElasticMaterial(pfem) {}
    
    //! calculate stress at material point
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! calculate tangent stiffness at material point
    tens4ds Tangent(FEMaterialPoint& pt) override;
    
    //! calculate strain energy density at material point
    double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! data initialization
    bool Validate() override;
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};
