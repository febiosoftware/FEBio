#include "stdafx.h"
#include "FEFacetSet.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FEFacetSet::FEFacetSet(FEMesh* mesh) : m_mesh(mesh)
{
}

//-----------------------------------------------------------------------------
void FEFacetSet::Create(int n)
{
	m_Face.resize(n);
}

//-----------------------------------------------------------------------------
FEFacetSet::FACET& FEFacetSet::Face(int i)
{
	return m_Face[i];
}

//-----------------------------------------------------------------------------
const FEFacetSet::FACET& FEFacetSet::Face(int i) const
{
	return m_Face[i];
}

//-----------------------------------------------------------------------------
void FEFacetSet::Add(FEFacetSet* pf)
{
	m_Face.insert(m_Face.end(), pf->m_Face.begin(), pf->m_Face.end());
}

//-----------------------------------------------------------------------------
FENodeSet FEFacetSet::GetNodeSet()
{
	FEMesh* pm = m_mesh;
	FENodeSet set(pm);

	vector<int> tag(pm->Nodes(), 0);
	for (int i = 0; i<Faces(); ++i)
	{
		FACET& el = m_Face[i];
		int ne = el.ntype;
		for (int j = 0; j<ne; ++j)
		{
			if (tag[el.node[j]] == 0)
			{
				set.add(el.node[j]);
				tag[el.node[j]] = 1;
			}
		}
	}
	return set;
}

//-----------------------------------------------------------------------------
void FEFacetSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_Face;
	}
	else
	{
		ar >> m_Face;
	}
}
