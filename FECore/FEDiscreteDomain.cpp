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
	m_Node.reserve(n);
	for (i=0; i<N; ++i) if (tag[i] >= 0) m_Node.push_back(i);
	assert(m_Node.size() == n);
	return true;
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
		for (int j=0; j<nint; ++j) el.m_State[j]->ShallowCopy(dmp, bsave);
	}
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
			ar << el.m_nrigid;
			ar << el.m_nID;
			ar << el.m_node;

			for (int j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
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
			ar >> el.m_nrigid;
			ar >> el.m_nID;
			ar >> el.m_node;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(fem.GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}
	}
}
