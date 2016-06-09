//
//  FEFiberMaterialPoint.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEFiberMaterialPoint.h"

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPoint* FEFiberMaterialPoint::Copy()
{
	FEFiberMaterialPoint* pt = new FEFiberMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEFiberMaterialPoint::Init()
{
	m_n0 = vec3d(1,0,0);

	// don't forget to intialize the base class data
	FEMaterialPoint::Init();
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEFiberMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
    
	if (ar.IsSaving())
	{
		ar << m_n0;
	}
	else
	{
		ar >> m_n0;
	}
}
