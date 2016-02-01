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
void FEEdge::create(int n)
{ 
	m_Elem.resize(n); 
}

//-----------------------------------------------------------------------------
//! serialization
void FEEdge::Serialize(DumpStream& ar)
{
}
