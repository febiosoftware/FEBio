#include "stdafx.h"
#include "FEPrescribedActiveContractionIsotropic.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPrescribedActiveContractionIsotropic, FEElasticMaterial)
	ADD_PARAMETER(m_T0 , "T0"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPrescribedActiveContractionIsotropic::FEPrescribedActiveContractionIsotropic(FEModel* pfem) : FEElasticMaterial(pfem)
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
