//
//  FEActiveConstantSupply.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 12/9/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEActiveConstantSupply.h"
#include "FEBioMech/FEElasticMaterial.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEActiveConstantSupply, FEActiveMomentumSupply)
    ADD_PARAMETER(m_asupp, FE_PARAM_DOUBLE, "supply");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEActiveConstantSupply::FEActiveConstantSupply(FEModel* pfem) : FEActiveMomentumSupply(pfem)
{
    m_asupp = 0;
}

//-----------------------------------------------------------------------------
//! Initialization.
void FEActiveConstantSupply::Init()
{
    FEMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Active momentum supply vector.
//! The momentum supply is oriented along the first material axis
vec3d FEActiveConstantSupply::ActiveSupply(FEMaterialPoint& mp)
{

    FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	
    // active momentum supply vector direction
    vec3d V(et.m_Q[0][0],et.m_Q[1][0], et.m_Q[2][0]);
    
    mat3d F = et.m_F;
    vec3d pw = (F*V)*m_asupp;

    return pw;
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
vec3d FEActiveConstantSupply::Tangent_ActiveSupply_Strain(FEMaterialPoint &mp)
{
    return vec3d(0,0,0);
}
