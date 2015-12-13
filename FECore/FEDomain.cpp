#include "stdafx.h"
#include "FEDomain.h"
#include "FEMaterial.h"
#include "FEDataExport.h"
#include "FEMesh.h"
#include <string.h>

//-----------------------------------------------------------------------------
FEDomain::FEDomain(int nclass, FEMesh* pm) : FECoreBase(FEDOMAIN_ID), m_pMesh(pm), m_nclass(nclass)
{
	m_szname[0] = 0;
}

//-----------------------------------------------------------------------------
FEDomain::~FEDomain()
{
	// delete all data export classes
	if (m_Data.empty() == false)
	{
		size_t ND = m_Data.size();
		for (size_t i=0; i<ND; ++i) delete m_Data[i];
		m_Data.clear();
	}
}

//-----------------------------------------------------------------------------
void FEDomain::AddDataExport(FEDataExport* pd)
{
	if (pd) m_Data.push_back(pd);
}

//-----------------------------------------------------------------------------
void FEDomain::SetName(const char* szname)
{
	if (szname == 0) return;
	int l = strlen(szname);
	if (l>=MAX_DOMAIN_NAME) l = MAX_DOMAIN_NAME - 1;
	if (l>0) strncpy(m_szname, szname, l);
	m_szname[l] = 0;
}

//-----------------------------------------------------------------------------
const char* FEDomain::GetName()
{
	return m_szname;
}

//-----------------------------------------------------------------------------
FEElement* FEDomain::FindElementFromID(int nid)
{
	for (int i=0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		if (el.m_nID == nid) return &el;
	}

	return 0;
}

//-----------------------------------------------------------------------------
void FEDomain::InitMaterialPointData()
{
	FEMaterial* pmat = GetMaterial();

	for (int i=0; i<Elements(); ++i)
	{
		FEElement& el = ElementRef(i);
		for (int k=0; k<el.GaussPoints(); ++k) el.SetMaterialPointData(pmat->CreateMaterialPointData(), k);
	}
}

//-----------------------------------------------------------------------------
void FEDomain::SetMatID(int mid)
{
	for (int i=0; i<Elements(); ++i) ElementRef(i).SetMatID(mid);
}

//-----------------------------------------------------------------------------
//! return a specific node
FENode& FEDomain::Node(int i)
{ 
	return m_pMesh->Node(m_Node[i]); 
}

//-----------------------------------------------------------------------------
bool FEDomain::Initialize(FEModel& fem)
{
	// make sure that there are elements in this domain
	if (Elements() == 0) return false;

	// get the mesh to which this domain belongs
	FEMesh& mesh = *GetMesh();

	// This array is used to keep tags on each node
	int NN = mesh.Nodes();
	vector<int> tag; tag.assign(NN, -1);

	// let's find all nodes the domain needs
	int nn = 0;
	int NE = Elements();
	for (int i=0; i<NE; ++i)
	{
		FEElement& el = ElementRef(i);
		int ne = el.Nodes();
		for (int j=0; j<ne; ++j)
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
	m_Node.assign(nn, -1);

	// fill the node index table
	for (int i=0; i<NN; ++i)
	{
		if (tag[i] >= 0)
		{
			m_Node[tag[i]] = i;
		}
	}

#ifdef _DEBUG
	// make sure all nodes are assigned a local index
	for (int i=0; i<nn; ++i)
	{
		assert(m_Node[i] >= 0);
	}
#endif

	return true;
}
