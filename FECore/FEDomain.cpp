#include "stdafx.h"
#include "FEDomain.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
FEDomain::FEDomain(int nclass, FEModel* fem) : FEMeshPartition(nclass, fem)
{

}

//-----------------------------------------------------------------------------
void FEDomain::SetMaterial(FEMaterial* pm)
{
	assert(pm);
	if (pm) pm->AddDomain(this);
}

//-----------------------------------------------------------------------------
void FEDomain::SetMatID(int mid)
{
	ForEachElement([=](FEElement& el) { el.SetMatID(mid); });
}

//-----------------------------------------------------------------------------
// This routine allocates the material point data for the element's integration points.
// Currently, this has to be called after the elements have been assigned a type (since this
// determines how many integration point an element gets). 
void FEDomain::CreateMaterialPointData()
{
	FEMaterial* pmat = GetMaterial();
	if (pmat) ForEachElement([=](FEElement& el) {
		for (int k = 0; k<el.GaussPoints(); ++k) el.SetMaterialPointData(pmat->CreateMaterialPointData(), k);
	});
}

//-----------------------------------------------------------------------------
// serialization
void FEDomain::Serialize(DumpStream& ar)
{
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			ar << m_Node;

			int NEL = Elements();
			for (int i = 0; i < NEL; ++i)
			{
				FEElement& el = ElementRef(i);
				el.Serialize(ar);
				int nint = el.GaussPoints();
				for (int j = 0; j < nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
		else
		{
			ar >> m_Node;

			FEMaterial* pmat = GetMaterial();
			int NEL = Elements();
			for (int i = 0; i < NEL; ++i)
			{
				FEElement& el = ElementRef(i);
				el.Serialize(ar);
				int nint = el.GaussPoints();
				for (int j = 0; j < nint; ++j)
				{
					el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
	else
	{
		int NEL = Elements();
		for (int i = 0; i < NEL; ++i)
		{
			FEElement& el = ElementRef(i);
			el.Serialize(ar);
			int nint = el.GaussPoints();
			for (int j = 0; j < nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
		}
	}
}
