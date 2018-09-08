#include "stdafx.h"
#include "FEEdge.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FEEdge::FEEdge(FEMesh* pm) : FEDomain(FE_DOMAIN_EDGE, pm)
{
}

//-----------------------------------------------------------------------------
FEEdge::~FEEdge()
{
}

//-----------------------------------------------------------------------------
FENodeSet FEEdge::GetNodeSet()
{
	FEMesh* pm = GetMesh();
	FENodeSet set(pm);

	vector<int> tag(pm->Nodes(), 0);
	for (int i = 0; i<Elements(); ++i)
	{
		FELineElement& el = Element(i);
		int ne = el.Nodes();
		for (int j = 0; j<ne; ++j)
		{
			if (tag[el.m_node[j]] == 0)
			{
				set.add(el.m_node[j]);
				tag[el.m_node[j]] = 1;
			}
		}
	}
	return set;
}

//-----------------------------------------------------------------------------
bool FEEdge::Init()
{
	// make sure that there is an edge defined
	if (Elements() == 0) return false;

	// get the mesh to which this edge belongs
	FEMesh& mesh = *GetMesh();

	// This array is used to keep tags on each node
	vector<int> tag; tag.assign(mesh.Nodes(), -1);

	// let's find all nodes the edge needs
	int nn = 0;
	int ne = Elements();
	for (int i=0; i<ne; ++i)
	{
		FELineElement& el = Element(i);
		el.m_lid = i;

		for (int j=0; j<el.Nodes(); ++j)
		{
			// get the global node number
			int m = el.m_node[j];
		
			// create a local node number
			if (tag[m] == -1) tag[m] = nn++;

			// set the local node number
			el.m_lnode[j] = tag[m];
		}
	}

	// allocate node index table
	m_Node.resize(nn);

	// fill the node index table
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		if (tag[i] >= 0)
		{
			m_Node[tag[i]] = i;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEEdge::Create(int nelems, int elemType)
{ 
	m_Elem.resize(nelems); 
	for (int i = 0; i < nelems; ++i)
	{
		FELineElement& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetDomain(this);
	}

	if (elemType != -1)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}
