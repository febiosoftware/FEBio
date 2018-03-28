#include "stdafx.h"
#include "FEDiscreteDomain.h"
#include "FEMesh.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
void FEDiscreteDomain::Create(int nelems, int elemType)
{ 
	m_Elem.resize(nelems); 
	if (elemType != -1)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
bool FEDiscreteDomain::Init()
{
	if (FEDomain::Init() == false) return false;

	FEMaterial* pmat = GetMaterial();
	if (pmat) SetMatID(pmat->GetID());

	return true;
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::AddElement(int eid, int n[2])
{
	FEDiscreteElement el;
	el.SetType(FE_DISCRETE);
	el.m_node[0] = n[0];
	el.m_node[1] = n[1];
	el.SetID(eid);
	m_Elem.push_back(el);
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::Serialize(DumpStream& ar)
{
	if (ar.IsShallow())
	{
		int NEL = (int) m_Elem.size();
		for (int i=0; i<NEL; ++i)
		{
			FEDiscreteElement& el = m_Elem[i];
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			ar << m_Node;
			int nel = (int) m_Elem.size();
			for (int i=0; i<nel; ++i)
			{
				FEDiscreteElement& el = m_Elem[i];
				int nmat = el.GetMatID();
				ar << (int) el.Type();
			
				ar << nmat;
				ar << el.GetID();
				ar << el.m_node;

				for (int j=0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
		else
		{
			FEMaterial* pmat = GetMaterial();
			assert(pmat);

			FEModel& fem = ar.GetFEModel();
			ar >> m_Node;
			int n, mat, nid;
			for (size_t i=0; i<m_Elem.size(); ++i)
			{
				FEDiscreteElement& el = m_Elem[i];
				ar >> n;

				el.SetType(n);

				ar >> mat; el.SetMatID(mat);
				ar >> nid; el.SetID(nid);
				ar >> el.m_node;

				for (int j=0; j<el.GaussPoints(); ++j)
				{
					el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
}
