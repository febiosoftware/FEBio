#include "stdafx.h"
#include "FEDiscreteDomain.h"
#include "FEMesh.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
bool FEDiscreteDomain::Initialize(FEModel &fem)
{
	int i, j;
	FEMesh& m = *m_pMesh;
	int N = m.Nodes();
	vector<int> tag; tag.assign(N, -1);

	int NE = Elements();
	int n = 0;
	for (i=0; i<NE; ++i)
	{
		FEDiscreteElement& e = m_Elem[i];
		int ne = e.Nodes();
		for (j=0; j<ne; ++j)
		{
			int nj = e.m_node[j];
			if (tag[nj] == -1) tag[nj] = n++;
		}
	}
	m_Node.resize(n);
	n = 0;
	for (i=0; i<NE; ++i)
	{
		FEDiscreteElement& e = m_Elem[i];
		int ne = e.Nodes();
		for (j=0; j<ne; ++j)
		{
			int nj = e.m_node[j];
			if (tag[nj] != -1)
			{
				m_Node[n++] = nj;
				tag[nj] = -1;
			}
		}
	}

	FEMaterial* pmat = GetMaterial();
	if (pmat)
	{
		int mid = pmat->GetID();
		for (int i=0; i<NE; ++i)
		{
			Element(i).SetMatID(mid);
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[DOF_X] = DOF_ACTIVE;
				node.m_ID[DOF_Y] = DOF_ACTIVE;
				node.m_ID[DOF_Z] = DOF_ACTIVE;
			}
		}
	}
}

//-----------------------------------------------------------------------------
FENode& FEDiscreteDomain::Node(int i) 
{
	return m_pMesh->Node(m_Node[i]); 
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::ShallowCopy(DumpStream& dmp, bool bsave)
{
	int NEL = (int) m_Elem.size();
	for (int i=0; i<NEL; ++i)
	{
		FEDiscreteElement& el = m_Elem[i];
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->ShallowCopy(dmp, bsave);
	}
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::AddElement(int eid, int n[2])
{
	FEDiscreteElement el;
	el.SetType(FE_DISCRETE);
	el.m_node[0] = n[0];
	el.m_node[1] = n[1];
	el.m_nID = eid;
	m_Elem.push_back(el);
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::Serialize(DumpFile& ar)
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
			ar << el.m_nID;
			ar << el.m_node;

			for (int j=0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
		}
	}
	else
	{
		FEModel& fem = *ar.GetFEModel();
		ar >> m_Node;
		int n, mat;
		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FEDiscreteElement& el = m_Elem[i];
			ar >> n;

			el.SetType(n);

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nID;
			ar >> el.m_node;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(fem.GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
	}
}
