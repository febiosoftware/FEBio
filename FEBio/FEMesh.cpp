// FEMesh.cpp: implementation of the FEMesh class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMesh.h"
#include "FEException.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEMesh::FEMesh()
{
	m_Domain.resize(3);
	m_Domain[0] = new FESolidDomain(this);
	m_Domain[1] = new FEShellDomain(this);
	m_Domain[2] = new FETrussDomain(this);
}

FEMesh::~FEMesh()
{

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
	m_Domain[0] = new FESolidDomain(this); *m_Domain[0] = *m.m_Domain[0];
	m_Domain[1] = new FESolidDomain(this); *m_Domain[1] = *m.m_Domain[1];
	m_Domain[2] = new FESolidDomain(this); *m_Domain[2] = *m.m_Domain[2];
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
	m_Domain[0] = new FESolidDomain(this); *m_Domain[0] = *m.m_Domain[0];
	m_Domain[1] = new FESolidDomain(this); *m_Domain[1] = *m.m_Domain[1];
	m_Domain[2] = new FESolidDomain(this); *m_Domain[2] = *m.m_Domain[2];

	return (*this);
}

//////////////////////////////////////////////////////////////////////
// FUNCTION: FEMesh::Create
//  Allocates storage for mesh data.
//

void FEMesh::Create(int nodes, int elems, int shells, int ntruss)
{
	if (nodes >0) m_Node.resize (nodes);

	if (elems >0) m_Domain[0]->create(elems);
	if (shells>0) m_Domain[1]->create(shells);
	if (ntruss>0) m_Domain[2]->create(ntruss);
}

int FEMesh::Elements()
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i) N += m_Domain[i]->Elements();
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
	vector<int> val(N);
	val.zero();

	// count the nodal valences
	FESolidDomain& bd = SolidDomain();
	for (i=0; i<bd.Elements(); ++i)
	{
		FEElement& el = bd.Element(i);
		n = el.Nodes();
		for (j=0; j<n; ++j) ++val[el.m_node[j]];
	}

	FEShellDomain& sd = ShellDomain();
	for (i=0; i<sd.Elements(); ++i)
	{
		FEElement& el = sd.Element(i);
		n = el.Nodes();
		for (j=0; j<n; ++j) ++val[el.m_node[j]];
	}

	FETrussDomain& td = TrussDomain();
	for (i=0; i<td.Elements(); ++i)
	{
		FEElement& el = td.Element(i);
		n = el.Nodes();
		for (j=0; j<n; ++j) ++val[el.m_node[j]];
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

		node.m_pt = 0;

		node.m_T = 0;
	}

	// update the mesh
	UpdateBox();

	// reset domain data
	for (int n=0; n<(int) m_Domain.size(); ++n) m_Domain[n]->Reset();
}

//-----------------------------------------------------------------------------
//! Update the mesh data. This function calculates the initial directors
//! for the shell elements.

// NOTE: This function needs to be called after the rigid bodies have been
// initialized

bool FEMesh::Init()
{
	int i, j, n, m0, m1, m2;
	int* en;
	vec3d a, b, c;

	vec3d* r0;

	int ninverted = 0;

	// check all solid elements to see if they are not initially inverted
	FESolidDomain& bd = SolidDomain();
	for (i=0; i<bd.Elements(); ++i)
	{
		FESolidElement& el = bd.Element(i);

		try
		{
			if (!el.IsRigid()) bd.UnpackElement(el);
		}
		catch (NegativeJacobian e)
		{
			fprintf(stderr, "**************************** E R R O R ****************************\n");
			fprintf(stderr, "Negative jacobian detected at integration point %d of element %d\n", e.m_ng, e.m_iel);
			fprintf(stderr, "Did you use the right node numbering?\n");
			if (e.m_pel)
			{
				FEElement &el = *e.m_pel;
				fprintf(stderr, "Nodes:");
				for (int n=0; n<el.Nodes(); ++n)
				{
					fprintf(stderr, "%d", el.m_node[n]+1);
					if (n+1 != el.Nodes()) fprintf(stderr, ","); else fprintf(stderr, "\n");
				}
			}
			fprintf(stderr, "*******************************************************************\n\n");
			++ninverted;
		}
	}


	// zero initial directors for shell nodes
	for (i=0; i<Nodes(); ++i) Node(i).m_D0 = vec3d(0,0,0);

	// initialize shell data
	FEShellDomain& sd = ShellDomain();
	for (i=0; i<sd.Elements(); ++i)
	{
		FEShellElement& el = sd.Element(i);
		sd.UnpackElement(el, 0);

		r0 = el.r0();

		n = el.Nodes();
		en = el.m_node;
		for (j=0; j<n; ++j)
		{

			m0 = j;
			m1 = (j+1)%n;
			m2 = (j==0? n-1: j-1);

			a = r0[m0];
			b = r0[m1];
			c = r0[m2];

			Node(en[m0]).m_D0 += (b-a)^(c-a);
		}
	}

	// make sure we start with unit directors
	for (i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		node.m_D0.unit();
		node.m_Dt = node.m_D0;
	}

	// check the connectivity of the shells
	for (i=0; i<sd.Elements(); ++i)
	{
		FEShellElement& el = sd.Element(i);

		try
		{
			if (!el.IsRigid()) sd.UnpackElement(el);
		}
		catch (NegativeJacobian e)
		{
			fprintf(stderr, "**************************** E R R O R ****************************\n");
			fprintf(stderr, "Negative jacobian detected at integration point %d of element %d\n", e.m_ng, e.m_iel);
			fprintf(stderr, "Did you use the right node numbering?\n");
			if (e.m_pel)
			{
				FEElement &el = *e.m_pel;
				fprintf(stderr, "Nodes:");
				for (int n=0; n<el.Nodes(); ++n)
				{
					fprintf(stderr, "%d", el.m_node[n]+1);
					if (n+1 != el.Nodes()) fprintf(stderr, ","); else fprintf(stderr, "\n");
				}
			}
			fprintf(stderr, "*******************************************************************\n\n");
			++ninverted;
		}
	}

	if (ninverted != 0)
	{
		fprintf(stderr, "**************************** E R R O R ****************************\n");
		fprintf(stderr, " FEBio found %d initially inverted elements.\n", ninverted);
		fprintf(stderr, " Run will be aborted.\n");
		fprintf(stderr, "*******************************************************************\n\n");
		return false;
	}

	// next if a node does not belong to a shell
	// we turn of the rotational degrees of freedom
	vector<int> tag(Nodes());
	tag.zero();
	for (i=0; i<sd.Elements(); ++i)
	{
		FEShellElement& el = sd.Element(i);
		if (!el.IsRigid()) sd.UnpackElement(el);

		n = el.Nodes();
		en = el.m_node;
		for (j=0; j<n; ++j) tag[en[j]] = 1;
	}

	for (i=0; i<Nodes(); ++i) 
	{
		FENode& node = Node(i);
		if (tag[i] == 0)
		{
			node.m_ID[3] = -1;
			node.m_ID[4] = -1;
			node.m_ID[5] = -1;
		}
	}

	for (i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		for (j=0; j<MAX_NDOFS; ++j)	node.m_BC[j] = node.m_ID[j];
	}

	// reset data
	Reset();

	return true;
}

//-----------------------------------------------------------------------------
//! This function calculates the volume of an element
//! using the element's integration points. It also assumes the element 
//! is unpacked.

double FEMesh::ElementVolume(FEElement& el)
{
	// determine the type of the element
	double V = 0;

	FESolidElement* ph = dynamic_cast<FESolidElement*>(&el);
	FESolidDomain& bd = SolidDomain();
	if (ph) 
	{
		bd.UnpackElement(*ph);
		int nint = ph->GaussPoints();
		double *w = ph->GaussWeights();
		for (int n=0; n<nint; ++n) V += ph->detJ0(n)*w[n];
	}

	FEShellDomain& sd = ShellDomain();
	FEShellElement* ps = dynamic_cast<FEShellElement*>(&el);
	if (ps) 
	{
		sd.UnpackElement(*ps);
		int nint = ps->GaussPoints();
		double *w = ps->GaussWeights();
		for (int n=0; n<nint; ++n) V += ps->detJ0(n)*w[n];
	}

	FESurfaceElement* pf = dynamic_cast<FESurfaceElement*>(&el);
	if (pf) V = 0;

	return V;
}

//-----------------------------------------------------------------------------
//! Assigns a material ID to the entire mesh

void FEMesh::SetMatID(int n)
{
	int i;
	FESolidDomain& bd = SolidDomain();
	int N = bd.Elements();
	for (i=0; i<N; ++i) bd.Element(i).SetMatID(n);

	FEShellDomain& sd = ShellDomain();
	N = sd.Elements();
	for (i=0; i<N; ++i) sd.Element(i).SetMatID(n);
}

//-----------------------------------------------------------------------------
//! Find a nodeset by ID

FENodeSet* FEMesh::FindNodeSet(int nid)
{
	for (int i=0; i<m_NodeSet.size(); ++i) if (m_NodeSet[i].GetID() == nid) return &m_NodeSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a nodeset by name

FENodeSet* FEMesh::FindNodeSet(const char* szname)
{
	for (int i=0; i<m_NodeSet.size(); ++i) if (strcmp(m_NodeSet[i].GetName(), szname)) return &m_NodeSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! This function is used by the restart feature and reads or writes
//! the mesh data to or from the binary archive
//! \param[in] ar the archive to which the data is serialized
//! \sa FEM::Serialize()
//! \sa Archive
//! \todo serialize nodesets

void FEMesh::Serialize(Archive& ar)
{
	if (ar.IsSaving())
	{
		int i;

		FESolidDomain& bd = SolidDomain();
		FEShellDomain& sd = ShellDomain();
		FETrussDomain& td = TrussDomain();

		// write mesh item counts
		int nn   = Nodes();
		int nbel = bd.Elements();
		int nsel = sd.Elements();
		int ntel = td.Elements();
		ar << nn << nbel << nsel << ntel;

		// write nodal data
		for (i=0; i<nn; ++i) ar.write(&Node(i), sizeof(FENode), 1);

		// write solid element data
		int nmat;
		for (i=0; i<nbel; ++i)
		{
			FESolidElement& el = bd.Element(i);
			nmat = el.GetMatID();
			ar << el.Type();
			
			ar << nmat;
			ar << el.m_nrigid;
			ar << el.m_nID;
			ar << el.m_node;

			ar << el.m_eJ;
			ar << el.m_ep;
			ar << el.m_Lk;
		}

		// write shell element data
		for (i=0; i<nsel; ++i)
		{
			FEShellElement& el = sd.Element(i);
			ar << el.Type();

			ar << el.m_eJ;
			ar << el.m_ep;

			ar << el.GetMatID();
			ar << el.m_nrigid;
			ar << el.m_nID;
			ar << el.m_node;

			ar << el.m_h0;
			ar << el.m_Lk;
		}

		// write bounding box data
		FE_BOUNDING_BOX box = GetBoundingBox();
		ar << box.r0 << box.r1;
	}
	else
	{
		int i, n, mat;

		// read mesh item counts
		int nn, nbel, nsel, ntel;
		ar >> nn >> nbel >> nsel >> ntel;

		// allocate storage for mesh data
		Create(nn, nbel, nsel, ntel);

		// read nodal data
		for (i=0; i<nn; ++i) ar.read(&Node(i), sizeof(FENode), 1);

		// read solid element data
		FESolidDomain& bd = SolidDomain();
		for (i=0; i<nbel; ++i)
		{
			FESolidElement& el = bd.Element(i);
			ar >> n;

			el.SetType(n);

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nrigid;
			ar >> el.m_nID;
			ar >> el.m_node;

			ar >> el.m_eJ;
			ar >> el.m_ep;
			ar >> el.m_Lk;
		}

		// read shell element data
		FEShellDomain& sd = ShellDomain();
		for (i=0; i<nsel; ++i)
		{
			FEShellElement& el = sd.Element(i);
			ar >> n;

			el.SetType(n);

			ar >> el.m_eJ;
			ar >> el.m_ep;

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nrigid;
			ar >> el.m_nID;
			ar >> el.m_node;

			ar >> el.m_h0;
			ar >> el.m_Lk;
		}

		// read bounding box data
		FE_BOUNDING_BOX& box = GetBoundingBox();
		ar >> box.r0 >> box.r1;
	}
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
	case FE_TET: return 4;
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
	int* en = el.m_node;
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

	FESolidDomain& bd = SolidDomain();
	FEShellDomain& sd = ShellDomain();
	FETrussDomain& td = TrussDomain();

	// search solid elements
	pe = bd.FindElementFromID(nid);

	// search shell elements
	if (pe == 0) sd.FindElementFromID(nid);

	// search truss elements
	if (pe == 0) td.FindElementFromID(nid);

	return pe;
}
