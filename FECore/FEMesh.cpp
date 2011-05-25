// FEMesh.cpp: implementation of the FEMesh class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMesh.h"
#include "FEException.h"
#include "FEDiscreteDomain.h"
#include "FETrussDomain.h"
#include "FEShellDomain.h"
#include "FESolidDomain.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEMesh::FEMesh()
{
}

FEMesh::~FEMesh()
{
	for (size_t i=0; i<m_NodeSet.size(); ++i) delete m_NodeSet[i];
	m_NodeSet.clear();
	ClearDomains();
}

void FEMesh::ClearDomains()
{
	for (int i=0; i<(int) m_Domain.size(); ++i) delete m_Domain[i];
	m_Domain.clear();
}

FEMesh::FEMesh(FEMesh& m)
{
	// copy nodal data
	m_Node = m.m_Node;

	// clear the domains
	ClearDomains();

	// copy domains
	for (int i=0; i<m.Domains(); ++i)
	{
		FEDomain* pd = m.Domain(i).Clone();
		pd->SetMesh(this);
		m_Domain.push_back(pd);
	}
}

//////////////////////////////////////////////////////////////////////
// FUNCTION: FEMesh::operator = 
//  Assignment operator.
//

FEMesh& FEMesh::operator =(FEMesh& m)
{
	// copy nodal data
	m_Node = m.m_Node;

	// clear the domains
	ClearDomains();

	// copy domains
	for (int i=0; i<m.Domains(); ++i)
	{
		FEDomain* pd = m.Domain(i).Clone();
		pd->SetMesh(this);
		m_Domain.push_back(pd);
	}

	return (*this);
}

//////////////////////////////////////////////////////////////////////
// FUNCTION: FEMesh::Create
//  Allocates storage for mesh data.
//

void FEMesh::CreateNodes(int nodes)
{
	assert(nodes);
	m_Node.resize (nodes);
}

void FEMesh::AddNode(vec3d r)
{
	FENode node;
	node.m_r0 = r;
	node.m_rt = node.m_r0;

	// set rigid body id
	node.m_rid = -1;

	// open displacement dofs
	node.m_ID[0] = 0;
	node.m_ID[1] = 0;
	node.m_ID[2] = 0;

	// open rotational dofs
	node.m_ID[3] = 0;
	node.m_ID[4] = 0;
	node.m_ID[5] = 0;

	// open pressure dof
	node.m_ID[6] = 0;

	// close the rigid rotational dofs
	node.m_ID[7] = -1;
	node.m_ID[8] = -1;
	node.m_ID[9] = -1;

	// fix temperature dof
	node.m_ID[10] = -1;

	// open concentration dof
	node.m_ID[11] = 0;

	m_Node.push_back(node);
}

int FEMesh::Elements()
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i) N += m_Domain[i]->Elements();
	return N;
}

int FEMesh::SolidElements()
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i)
	{
		FESolidDomain* pd = dynamic_cast<FESolidDomain*>(m_Domain[i]);
		if (pd) N += pd->Elements();
	}
	return N;
}

int FEMesh::ShellElements()
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i)
	{
		FEShellDomain* pd = dynamic_cast<FEShellDomain*>(m_Domain[i]);
		if (pd) N += pd->Elements();
	}
	return N;
}

int FEMesh::TrussElements()
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i)
	{
		FETrussDomain* pd = dynamic_cast<FETrussDomain*>(m_Domain[i]);
		if (pd) N += pd->Elements();
	}
	return N;
}

int FEMesh::DiscreteElements()
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i)
	{
		FEDiscreteDomain* pd = dynamic_cast<FEDiscreteDomain*>(m_Domain[i]);
		if (pd) N += pd->Elements();
	}
	return N;
}


///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEMesh::UpdateBox
//  Updates the bounding box of the mesh (using current coordinates)
//

void FEMesh::UpdateBox()
{
	vec3d r0, r1;

	r0 = r1 = Node(0).m_rt;
	for (int i=1; i<Nodes(); ++i)
	{
		vec3d r = Node(i).m_rt;

		if (r.x < r0.x) r0.x = r.x;
		if (r.y < r0.y) r0.y = r.y;
		if (r.z < r0.z) r0.z = r.z;

		if (r.x > r1.x) r1.x = r.x;
		if (r.y > r1.y) r1.y = r.y;
		if (r.z > r1.z) r1.z = r.z;
	}

	m_box.r0 = r0;
	m_box.r1 = r1;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEMesh::RemoveIsolatedVertices
//  Counts the number of shell elements in the mesh
//

int FEMesh::RemoveIsolatedVertices()
{
	int i, j, k, N = Nodes(), n;

	// create a valence array
	vector<int> val; val.assign(N, 0);

	// count the nodal valences
	for (i=0; i<(int) m_Domain.size(); ++i)
	{
		FEDomain& d = Domain(i);
		for (j=0; j<d.Elements(); ++j)
		{
			FEElement& el = d.ElementRef(j);
			n = el.Nodes();
			for (k=0; k<n; ++k) ++val[el.m_node[k]];
		}
	}

	// see if there are any isolated nodes
	int ni = 0;
	for (i=0; i<N; ++i)
		if (val[i] == 0)
		{
			++ni;
			FENode& node = Node(i);
			for (k=0; k<MAX_NDOFS; ++k) node.m_ID[k] = -1;
		}

	return ni;
}

//-----------------------------------------------------------------------------
//! Reset the mesh data. Return nodes to their intial position, reset their 
//! attributes and zero all element stresses.

void FEMesh::Reset()
{
	// reset nodal data
	for (int i=0; i<Nodes(); ++i) 
	{
		FENode& node = Node(i);

		node.m_rp = node.m_rt = node.m_r0;
		node.m_vp = node.m_vt = node.m_v0;
		node.m_ap = node.m_at = vec3d(0,0,0);

		node.m_pt = node.m_p0;
		
		node.m_ct = node.m_cp = node.m_c0;
		
		node.m_T = 0;
	}

	// update the mesh
	UpdateBox();

	// reset domain data
	for (int n=0; n<(int) m_Domain.size(); ++n) m_Domain[n]->Reset();
}

//-----------------------------------------------------------------------------
//! This function calculates the volume of an element
//! using the element's integration points. It also assumes the element 
//! is unpacked.

double FEMesh::ElementVolume(FEElement& el)
{
	// determine the type of the element
	double V = 0;

	// get the domain from the domain ID of the element
	FEDomain* pd = m_Domain[el.m_gid];

	if (dynamic_cast<FESolidDomain*>(pd))
	{
		FESolidElement* ph = dynamic_cast<FESolidElement*>(&el);
		FESolidDomain& bd = dynamic_cast<FESolidDomain&>(*pd);
		pd->UnpackElement(*ph);
		int nint = ph->GaussPoints();
		double *w = ph->GaussWeights();
		for (int n=0; n<nint; ++n) V += bd.detJ0(*ph, n)*w[n];
	}

	if (dynamic_cast<FEShellDomain*>(pd))
	{
		FEShellElement* ps = dynamic_cast<FEShellElement*>(&el);
		FEShellDomain& sd = dynamic_cast<FEShellDomain&>(*pd);
		pd->UnpackElement(*ps);
		int nint = ps->GaussPoints();
		double *w = ps->GaussWeights();
		for (int n=0; n<nint; ++n) V += sd.detJ0(*ps, n)*w[n];
	}

	FESurfaceElement* pf = dynamic_cast<FESurfaceElement*>(&el);
	if (pf) V = 0;

	return V;
}

//-----------------------------------------------------------------------------
//! Find a nodeset by ID

FENodeSet* FEMesh::FindNodeSet(int nid)
{
	for (size_t i=0; i<m_NodeSet.size(); ++i) if (m_NodeSet[i]->GetID() == nid) return m_NodeSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a nodeset by name

FENodeSet* FEMesh::FindNodeSet(const char* szname)
{
	for (size_t i=0; i<m_NodeSet.size(); ++i) if (strcmp(m_NodeSet[i]->GetName(), szname)) return m_NodeSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
int FEMesh::Faces(FEElement& el)
{
	switch (el.Type())
	{
	case FE_HEX:
	case FE_RIHEX:
	case FE_UDGHEX: return 6;
	case FE_PENTA: return 5;
	case FE_TET: 
	case FE_TETG1: return 4;
	case FE_SHELL_QUAD:
	case FE_SHELL_TRI: return 1;
	}

	return 0;
}

//-----------------------------------------------------------------------------
//! This function returns the face connectivity from a certain element

int FEMesh::GetFace(FEElement& el, int n, int nf[4])
{
	int nn = -1;
	int* en = &el.m_node[0];
	switch (el.Type())
	{
	case FE_HEX:
	case FE_RIHEX:
	case FE_UDGHEX:
		nn = 4;
		switch (n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; break;
		case 3: nf[0] = en[0]; nf[1] = en[4]; nf[2] = en[7]; nf[3] = en[3]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; break;
		}
		break;
	case FE_PENTA:
		switch(n)
		{
		case 0: nn = 4; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; nf[3] = en[3]; break;
		case 1: nn = 4; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[5]; nf[3] = en[4]; break;
		case 2: nn = 4; nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[5]; nf[3] = en[2]; break;
		case 3: nn = 3; nf[0] = en[0]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[1]; break;
		case 4: nn = 3; nf[0] = en[3]; nf[1] = en[4]; nf[2] = en[5]; nf[3] = en[5]; break;
		}
		break;
	case FE_TET:
	case FE_TETG1:
		nn = 3;
		switch (n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = nf[3] = en[3]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = nf[3] = en[3]; break;
		case 2: nf[0] = en[0]; nf[1] = en[3]; nf[2] = nf[3] = en[2]; break;
		case 3: nf[0] = en[0]; nf[1] = en[2]; nf[2] = nf[3] = en[1]; break;
		}
		break;
	case FE_SHELL_QUAD:
		nn = 4;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3];
		break;
	case FE_SHELL_TRI:
		nn = 3;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2];
		break;
	}

	return nn;
}

//-----------------------------------------------------------------------------
//! Find an element from a given ID. return 0 if the element cannot be found.

FEElement* FEMesh::FindElementFromID(int nid)
{
	FEElement* pe = 0;

	for (int i=0; i<Domains(); ++i)
	{
		FEDomain& d = Domain(i);
		pe = d.FindElementFromID(nid);
		if (pe) return pe;
	}

	return pe;
}

//-----------------------------------------------------------------------------
// Find the element in which point y lies
FESolidElement* FEMesh::FindSolidElement(vec3d y, double r[3])
{
	int ND = (int) m_Domain.size();
	for (int i=0; i<ND; ++i)
	{
		FESolidDomain* pd = dynamic_cast<FESolidDomain*>(m_Domain[i]);
		if (pd)
		{
			FESolidElement* pe = pd->FindElement(y, r);
			if (pe) return pe;
		}
	}
	return 0;
}
