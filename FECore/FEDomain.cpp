#include "stdafx.h"
#include "FEDomain.h"
#include "FEMaterial.h"
#include "FEDataExport.h"
#include "FEMesh.h"
#include "DOFS.h"
#include "FEGlobalMatrix.h"
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
void FEDomain::SetDOF(vector<int>& dof)
{
	m_dof = dof;
}

//-----------------------------------------------------------------------------
// This is the default packing method. 
// It stores all the degrees of freedom for the first node in the order defined
// by the DOF array, then for the second node, and so on. 
void FEDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	int ndofs = m_dof.size();
	lm.resize(N*ndofs);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;
		for (int j=0; j<ndofs; ++j) lm[i*ndofs + j] = id[m_dof[j]];
	}
}

//-----------------------------------------------------------------------------
void FEDomain::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> elm;
	const int NE = Elements();
	for (int j=0; j<NE; ++j)
	{
		FEElement& el = ElementRef(j);
		UnpackLM(el, elm);
		M.build_add(elm);
	}
}

//-----------------------------------------------------------------------------
void FEDomain::Activate()
{
	// get the number of degrees of freedom for this domain.
	const int ndofs = (int)m_dof.size();

	// activate all the degrees of freedom of this domain
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.m_bexclude == false)
		{
			for (int j=0; j<ndofs; ++j) node.m_ID[m_dof[j]] = DOF_ACTIVE;
		}
	}
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
