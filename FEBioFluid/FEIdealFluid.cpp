//
//  FEIdealFluid.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/16/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEIdealFluid.h"
#include "FEFluid.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEIdealFluid, FEElasticFluid)
    ADD_PARAMETER2(m_k, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "k");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FEIdealFluid::FEIdealFluid(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_k = 1;
}

//-----------------------------------------------------------------------------
//! fluid pressure (relative to reference pressure when J=1)
double FEIdealFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
    return m_k*(1 - pt.m_J);
}

//-----------------------------------------------------------------------------
//! tangent of fluid pressure with respect to strain J
double FEIdealFluid::Tangent_Pressure_Strain(FEMaterialPoint &mp)
{
    return -m_k;
}
