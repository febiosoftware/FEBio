#include "stdafx.h"
#include "FEDiscreteMaterial.h"

//=============================================================================
FEMaterialPoint* FEDiscreteMaterialPoint::Copy()
{
	FEDiscreteMaterialPoint* pt = new FEDiscreteMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEDiscreteMaterialPoint::Serialize(DumpFile& ar)
{
	if (m_pNext) m_pNext->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEDiscreteMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (m_pNext) m_pNext->ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FEDiscreteMaterialPoint::Init(bool bflag)
{
	if (m_pNext) m_pNext->Init(bflag);
}

//=============================================================================
FEDiscreteMaterial::FEDiscreteMaterial(FEModel* pfem) : FEMaterial(pfem) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEDiscreteMaterial::CreateMaterialPointData() { return new FEDiscreteMaterialPoint; }
