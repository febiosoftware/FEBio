//
//  FEDamageMaterialPoint.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEDamageMaterialPoint.h"
#include <FECore/DumpStream.h>

FEMaterialPoint* FEDamageMaterialPoint::Copy()
{
    FEDamageMaterialPoint* pt = new FEDamageMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

void FEDamageMaterialPoint::Init()
{
	FEMaterialPoint::Init();

	// intialize data to zero
	m_Emax = 0;
	m_Etrial = 0;
	m_D = 0;
}

void FEDamageMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPoint::Update(timeInfo);

	m_Emax = max(m_Emax, m_Etrial);
}

void FEDamageMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	ar & m_Etrial & m_Emax & m_D;
}
