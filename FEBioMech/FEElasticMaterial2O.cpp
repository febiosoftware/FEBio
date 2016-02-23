#include "stdafx.h"
#include "FEElasticMaterial2O.h"
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
//! constructor
FEElasticMaterialPoint2O::FEElasticMaterialPoint2O()
{
}

//-----------------------------------------------------------------------------
//! initialization
void FEElasticMaterialPoint2O::init(bool bflag)
{
	if (bflag)
	{
	}
	else
	{
	}
}

//-----------------------------------------------------------------------------
//! create a shallow copy
FEMaterialPoint* FEElasticMaterialPoint2O::Copy()
{
	FEElasticMaterialPoint2O* pt = new FEElasticMaterialPoint2O(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! serialize material point data
void FEElasticMaterialPoint2O::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
	}
	else
	{
	}
}

//-----------------------------------------------------------------------------
FEElasticMaterial2O::FEElasticMaterial2O(FEModel* pfem) : FEElasticMaterial(pfem)
{
}
