#include "FEFiberMaterialPoint.h"

//! constructor
FEFiberMaterialPoint::FEFiberMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) 
{
	m_n0 = vec3d(1,0,0);
}

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPoint* FEFiberMaterialPoint::Copy()
{
	FEFiberMaterialPoint* pt = new FEFiberMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	pt->m_n0 = m_n0;
	return pt;
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
