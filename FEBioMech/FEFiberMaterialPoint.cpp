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
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEFiberMaterialPoint::Init(bool bflag)
{
	if (bflag)
	{
        m_n0 = vec3d(1,0,0);
	}
    
	// don't forget to intialize the nested data
	if (m_pt) m_pt->Init(bflag);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEFiberMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (m_pt) m_pt->ShallowCopy(dmp, bsave);
    
	if (bsave)
	{
		dmp << m_n0;
	}
	else
	{
		dmp >> m_n0;
	}
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEFiberMaterialPoint::Serialize(DumpFile& ar)
{
	if (m_pt) m_pt->Serialize(ar);
    
	if (ar.IsSaving())
	{
		ar << m_n0;
	}
	else
	{
		ar >> m_n0;
	}
}
