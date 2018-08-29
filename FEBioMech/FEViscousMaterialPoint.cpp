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

void FEViscousMaterialPoint::Init()
{
	// initialize data to identity
    m_Fp = mat3dd(1);

	m_dt = 0.0;
    
    // don't forget to intialize the nested data
    FEMaterialPoint::Init();
}

void FEViscousMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEElasticMaterialPoint& pt = *this->ExtractData<FEElasticMaterialPoint>();
    m_Fp = pt.m_F;

	m_dt = timeInfo.timeIncrement;
    
    // don't forget to call base class
    FEMaterialPoint::Update(timeInfo);
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
    mat3d L = (mat3dd(1) - m_Fp*Fi)/m_dt;
    return L;
}

mat3ds FEViscousMaterialPoint::RateOfDeformation()
{
    mat3d L = VelocityGradient();
    return L.sym();
}
