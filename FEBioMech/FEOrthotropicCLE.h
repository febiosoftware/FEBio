//
//  FEOrthotropicCLE.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/25/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements an orthotropic conewise linear elastic (CLE) material (Curnier et al. 1995 J Elasticity).
class FEOrthotropicCLE :	public FEElasticMaterial
{
public:
    double	lp11, lp22, lp33;	// diagonal first lamé constants (tension)
    double	lm11, lm22, lm33;	// diagonal first lamé constants (compression)
    double	l12, l23, l31;	// off-diagonal first lamé constants
    double	mu1, mu2, mu3;	// shear moduli
    
public:
    FEOrthotropicCLE(FEModel* pfem) : FEElasticMaterial(pfem) {}
    
    //! calculate stress at material point
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! calculate tangent stiffness at material point
    tens4ds Tangent(FEMaterialPoint& pt);
    
    //! calculate strain energy density at material point
    double StrainEnergyDensity(FEMaterialPoint& pt);
    
    //! data initialization
    bool Validate();
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};
