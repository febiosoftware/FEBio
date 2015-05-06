//
//  FEPrescribedActiveContractionIsotropic.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/13/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEPrescribedActiveContractionIsotropic.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEPrescribedActiveContractionIsotropic, FEElasticMaterial)
ADD_PARAMETER(m_T0 , FE_PARAM_DOUBLE, "T0"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPrescribedActiveContractionIsotropic::FEPrescribedActiveContractionIsotropic(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
void FEPrescribedActiveContractionIsotropic::Init()
{
}

//-----------------------------------------------------------------------------
mat3ds FEPrescribedActiveContractionIsotropic::Stress(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

    double J = pt.m_J;
    mat3ds b = pt.LeftCauchyGreen();
    
    // evaluate the active stress
    mat3ds s = b*(m_T0/J);
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPrescribedActiveContractionIsotropic::Tangent(FEMaterialPoint &mp)
{
    tens4ds c;
    c.zero();
    
    return c;
}
