//
//  FEPrescribedActiveContractionTransIso.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/13/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEPrescribedActiveContractionTransIso.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEPrescribedActiveContractionTransIso, FEElasticMaterial)
ADD_PARAMETER(m_T0 , FE_PARAM_DOUBLE, "T0"   );
ADD_PARAMETER(m_thd, FE_PARAM_DOUBLE, "theta");
ADD_PARAMETER(m_phd, FE_PARAM_DOUBLE, "phi"  );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPrescribedActiveContractionTransIso::FEPrescribedActiveContractionTransIso(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_thd = 0;
    m_phd = 90;
}

//-----------------------------------------------------------------------------
void FEPrescribedActiveContractionTransIso::Init()
{
    // convert angles from degrees to radians
    double pi = 4*atan(1.0);
    double the = m_thd*pi/180.;
    double phi = m_phd*pi/180.;
    // fiber direction in local coordinate system (reference configuration)
    m_n0.x = cos(the)*sin(phi);
    m_n0.y = sin(the)*sin(phi);
    m_n0.z = cos(phi);
}

//-----------------------------------------------------------------------------
mat3ds FEPrescribedActiveContractionTransIso::Stress(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    
    // get the initial fiber direction
    vec3d n0, nt;
    
    // evaluate fiber direction in global coordinate system
    n0 = pt.m_Q*m_n0;
    
    // evaluate the deformed fiber direction
    nt = F*n0;
    nt.unit();
    mat3ds N = dyad(nt);
    mat3dd I(1);
    
    // evaluate the active stress
    mat3ds s = (I-N)*m_T0;
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPrescribedActiveContractionTransIso::Tangent(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    
    // get the initial fiber direction
    vec3d n0, nt;
    
    // evaluate fiber direction in global coordinate system
    n0 = pt.m_Q*m_n0;
    
    // evaluate the deformed fiber direction
    nt = F*n0;
    nt.unit();
    mat3ds N = dyad(nt);
    mat3dd I(1);
    
    tens4ds c = (dyad1s(I) - dyad4s(I)*2 + dyad1s(N)*2 - dyad1s(I, N)/2.0)*m_T0;
    
    return c;
}
