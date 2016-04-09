//
//  FEViscousMaterialPoint.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 4/7/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FEViscousMaterialPoint.h"
#include "FEElasticMaterial.h"
#include <FECore/DumpStream.h>

FEMaterialPoint* FEViscousMaterialPoint::Copy()
{
    FEViscousMaterialPoint* pt = new FEViscousMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

void FEViscousMaterialPoint::Init(bool bflag)
{
    if (bflag)
    {
        // initialize data to identity
        m_Fp = mat3dd(1);
    }
    else
    {
        FEElasticMaterialPoint& pt = *this->ExtractData<FEElasticMaterialPoint>();
        m_Fp = pt.m_F;
    }
    
    // don't forget to intialize the nested data
    if (m_pNext) m_pNext->Init(bflag);
}

void FEViscousMaterialPoint::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << m_Fp;
    }
    else
    {
        ar >> m_Fp;
    }
    FEMaterialPoint::Serialize(ar);
}

mat3d FEViscousMaterialPoint::VelocityGradient()
{
    FEElasticMaterialPoint& pt = *this->ExtractData<FEElasticMaterialPoint>();
    mat3d Fi = pt.m_F.inverse();
    mat3d L = (mat3dd(1) - m_Fp*Fi)/dt;
    return L;
}

mat3ds FEViscousMaterialPoint::RateOfDeformation()
{
    mat3d L = VelocityGradient();
    return L.sym();
}
