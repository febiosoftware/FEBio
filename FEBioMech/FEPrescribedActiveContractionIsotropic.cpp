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
    mat3dd I(1);
    
    // evaluate the active stress
    mat3ds s = I*m_T0;
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPrescribedActiveContractionIsotropic::Tangent(FEMaterialPoint &mp)
{
    mat3dd I(1);
    
    tens4ds c = (dyad1s(I) - dyad4s(I)*2)*m_T0;
    
    return c;
}
