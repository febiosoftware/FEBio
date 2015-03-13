//
//  FEPrescribedActiveContractionIsotropicUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/13/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEPrescribedActiveContractionIsotropicUC.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEPrescribedActiveContractionIsotropicUC, FEUncoupledMaterial)
ADD_PARAMETER(m_T0 , FE_PARAM_DOUBLE, "T0"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPrescribedActiveContractionIsotropicUC::FEPrescribedActiveContractionIsotropicUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
void FEPrescribedActiveContractionIsotropicUC::Init()
{
}

//-----------------------------------------------------------------------------
// Since the prescribed active contraction stress is not dependent on a strain
// energy density function, we don't return the deviatoric part of the
// stress.  Instead, we return the actual stress.
mat3ds FEPrescribedActiveContractionIsotropicUC::DevStress(FEMaterialPoint &mp)
{
    mat3dd I(1);
    
    // evaluate the active stress
    mat3ds s = I*m_T0;
    
    return s;
}

//-----------------------------------------------------------------------------
// Since the prescribed active contraction stress is not dependent on a strain
// energy density function, we don't return the deviatoric part of the
// tangent.  Instead, we return the actual tangent.
tens4ds FEPrescribedActiveContractionIsotropicUC::DevTangent(FEMaterialPoint &mp)
{
    mat3dd I(1);
    
    tens4ds c = (dyad1s(I) - dyad4s(I)*2)*m_T0;
    
    return c;
}
