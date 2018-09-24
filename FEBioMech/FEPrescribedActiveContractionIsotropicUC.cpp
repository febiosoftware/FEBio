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
	ADD_PARAMETER(m_T0 , "T0"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPrescribedActiveContractionIsotropicUC::FEPrescribedActiveContractionIsotropicUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
mat3ds FEPrescribedActiveContractionIsotropicUC::DevStress(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    double J = pt.m_J;
    mat3ds b = pt.LeftCauchyGreen();
    
    // evaluate the active stress
    mat3ds s = b*(m_T0/J);
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPrescribedActiveContractionIsotropicUC::DevTangent(FEMaterialPoint &mp)
{
    return tens4ds(0.0);
}
