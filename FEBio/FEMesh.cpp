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
}

FEMesh::~FEMesh()
{

}

FEMesh::FEMesh(FEMesh& m)
{
	// copy nodal data
	m_Node = m.m_Node;

	// copy element data
	m_Elem  = m.m_Elem;
	m_Shell = m.m_Shell;
}

//////////////////////////////////////////////////////////////////////
// FUNCTION: FEMesh::operator = 
//  Assignment operator.
//

FEMesh& FEMesh::operator =(FEMesh& m)
{
	// copy nodal data
	m_Node = m.m_Node;

	// copy element data
	m_Elem  = m.m_Elem;
	m_Shell = m.m_Shell;

	return (*this);
}

//////////////////////////////////////////////////////////////////////
// FUNCTION: FEMesh::Create
//  Allocates storage for mesh data.
//

void FEMesh::Create(int nodes, int elems, int shells)
{
	if (nodes >0) m_Node.create (nodes);
	if (elems >0) m_Elem.create (elems);
	if (shells>0) m_Shell.create(shells);
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEMesh::UpdateBox
//  Updates the bounding box of the mesh (using current coordinates)
//

void FEMesh::UpdateBox()
{
	vec3d r0, r1;

	for (int i=0; i<Nodes(); ++i)
	{
		vec3d r = Node(i).m_rt;

		if (r.x < r0.x) r0.x = r.x;
		if (r.y < r0.y) r0.y = r.y;
		if (r.z < r0.z) r0.z = r.z;

		if (r.x > r0.x) r1.x = r.x;
		if (r.y > r0.y) r1.y = r.y;
		if (r.z > r0.z) r1.z = r.z;
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
	for (i=0; i<SolidElements(); ++i)
	{
		FEElement& el = SolidElement(i);
		n = el.Nodes();
		for (j=0; j<n; ++j) ++val[el.m_node[j]];
	}

	for (i=0; i<ShellElements(); ++i)
	{
		FEElement& el = ShellElement(i);
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
//! Unpack the element. That is, copy element data in traits structure

void FEMesh::UnpackElement(FESolidElement& el, unsigned int nflag)
{
	int i, n;

	vec3d* rt = el.rt();
	vec3d* r0 = el.r0();
	vec3d* vt = el.vt();
	double* pt = el.pt();

	int N = el.Nodes();
	int* lm = el.LM();

	for (i=0; i<N; ++i)
	{
		n = el.m_node[i];
		FENode& node = Node(n);

		int* id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[0];
		lm[3*i+1] = id[1];
		lm[3*i+2] = id[2];

		// now the pressure dofs
		lm[3*N+i] = id[6];

		// rigid rotational dofs
		lm[4*N + 3*i  ] = id[7];
		lm[4*N + 3*i+1] = id[8];
		lm[4*N + 3*i+2] = id[9];

		// fill the rest with -1
		lm[7*N + 3*i  ] = -1;
		lm[7*N + 3*i+1] = -1;
		lm[7*N + 3*i+2] = -1;
	}

	// copy nodal data to element arrays
	for (i=0; i<N; ++i)
	{
		n = el.m_node[i];

		FENode& node = Node(n);

		// initial coordinates (= material coordinates)
		r0[i] = node.m_r0;

		// current coordinates (= spatial coordinates)
		rt[i] = node.m_rt;

		// current nodal pressures
		pt[i] = node.m_pt;

		// current nodal velocities
		vt[i] = node.m_vt;
	}

	// unpack the traits data
	el.UnpackTraitsData(nflag);
}

//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure
//! Note that for the shell elements the lm order is different compared
//! to the solid element ordering. This is because for shell elements the
//! nodes have six degrees of freedom each, where for solids they only
//! have 3 dofs.

void FEMesh::UnpackElement(FEShellElement& el, unsigned int nflag)
{
	int i, n;

	vec3d* rt = el.rt();
	vec3d* r0 = el.r0();
	vec3d* vt = el.vt();
	vec3d* D0 = el.D0();
	vec3d* Dt = el.Dt();
	double* pt = el.pt();

	int N = el.Nodes();
	int* lm = el.LM();

	for (i=0; i<N; ++i)
	{
		n = el.m_node[i];
		FENode& node = Node(n);

		int* id = node.m_ID;

		// first the displacement dofs
		lm[6*i  ] = id[0];
		lm[6*i+1] = id[1];
		lm[6*i+2] = id[2];

		// next the rotational dofs
		lm[6*i+3] = id[3];
		lm[6*i+4] = id[4];
		lm[6*i+5] = id[5];

		// now the pressure dofs
		lm[6*N+i] = id[6];

		// rigid rotational dofs
		lm[7*N + 3*i  ] = id[7];
		lm[7*N + 3*i+1] = id[8];
		lm[7*N + 3*i+2] = id[9];

	}

	// copy nodal data to element arrays
	for (i=0; i<N; ++i)
	{
		n = el.m_node[i];

		FENode& node = Node(n);

		// initial coordinates (= material coordinates)
		r0[i] = node.m_r0;

		// current coordinates (= spatial coordinates)
		rt[i] = node.m_rt;

		// current nodal pressures
		pt[i] = node.m_pt;

		// current nodal velocities
		vt[i] = node.m_vt;

		// intial director
		D0[i] = node.m_D0;

		// current director
		Dt[i] = node.m_Dt;
	}

	// unpack the traits data
	el.UnpackTraitsData(nflag);
}

//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure

void FEMesh::UnpackElement(FESurfaceElement& el, unsigned int nflag)
{
	int i, n;

	vec3d* rt = el.rt();
	vec3d* r0 = el.r0();
	vec3d* vt = el.vt();
	double* pt = el.pt();

	int N = el.Nodes();
	int* lm = el.LM();

	for (i=0; i<N; ++i)
	{
		n = el.m_node[i];
		FENode& node = Node(n);

		int* id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[0];
		lm[3*i+1] = id[1];
		lm[3*i+2] = id[2];

		// now the pressure dofs
		lm[3*N+i] = id[6];

		// rigid rotational dofs
		lm[4*N + 3*i  ] = id[7];
		lm[4*N + 3*i+1] = id[8];
		lm[4*N + 3*i+2] = id[9];

		// fill the rest with -1
		lm[7*N + 3*i  ] = -1;
		lm[7*N + 3*i+1] = -1;
		lm[7*N + 3*i+2] = -1;
	}

	// copy nodal data to element arrays
	for (i=0; i<N; ++i)
	{
		n = el.m_node[i];

		FENode& node = Node(n);

		// initial coordinates (= material coordinates)
		r0[i] = node.m_r0;

		// current coordinates (= spatial coordinates)
		rt[i] = node.m_rt;

		// current nodal pressures
		pt[i] = node.m_pt;

		// current nodal velocities
		vt[i] = node.m_vt;
	}

	// unpack the traits data
	el.UnpackTraitsData(nflag);
}

//-----------------------------------------------------------------------------
//! Reset the mesh data. Return nodes to their intial position, reset their 
//! attributes and zero all element stresses.

void FEMesh::Reset()
{
	int i;

	// reset nodal data
	for (i=0; i<Nodes(); ++i) 
	{
		FENode& node = Node(i);

		node.m_rp = node.m_rt = node.m_r0;
		node.m_vp = node.m_vt = node.m_v0;
		node.m_ap = node.m_at = vec3d(0,0,0);

		node.m_pt = 0;
	}

	// update the mesh
	UpdateBox();

	// initialize solid element data
	for (i=0; i<SolidElements(); ++i) 
	{
		FESolidElement& el = SolidElement(i);
		el.Init(true);
	}

	// initialize shell element data
	for (i=0; i<ShellElements(); ++i)
	{
		FEShellElement& el = ShellElement(i);
		el.Init(true);
	}
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
	for (i=0; i<SolidElements(); ++i)
	{
		FESolidElement& el = SolidElement(i);

		try
		{
			if (!el.IsRigid()) UnpackElement(el);
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
					fprintf(stderr, "%d", el.m_node[n]);
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
	for (i=0; i<ShellElements(); ++i)
	{
		FEShellElement& el = ShellElement(i);
		UnpackElement(el, 0);

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
	for (i=0; i<ShellElements(); ++i)
	{
		FEShellElement& el = ShellElement(i);

		try
		{
			if (!el.IsRigid()) UnpackElement(el);
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
					fprintf(stderr, "%d", el.m_node[n]);
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
	for (i=0; i<ShellElements(); ++i)
	{
		FEShellElement& el = ShellElement(i);
		if (!el.IsRigid()) UnpackElement(el);

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
	if (ph) 
	{
		UnpackElement(*ph);
		int nint = ph->GaussPoints();
		double *w = ph->GaussWeights();
		for (int n=0; n<nint; ++n) V += ph->detJ0(n)*w[n];
	}

	FEShellElement* ps = dynamic_cast<FEShellElement*>(&el);
	if (ps) 
	{
		UnpackElement(*ps);
		int nint = ps->GaussPoints();
		double *w = ps->GaussWeights();
		for (int n=0; n<nint; ++n) V += ps->detJ0(n)*w[n];
	}

	FESurfaceElement* pf = dynamic_cast<FESurfaceElement*>(&el);
	if (pf) V = 0;

	return V;
}

//-----------------------------------------------------------------------------
//! Calculates the exact volume of a hexahedral element

double FEMesh::HexVolume(FESolidElement& el)
{
	// let's make sure this is indeed a hex element
	assert(el.Type() == FE_HEX);

	// get the nodal coordinates
	vec3d* r = el.rt();
	double x1 = r[0].x, y1 = r[0].y, z1 = r[0].z;
	double x2 = r[1].x, y2 = r[1].y, z2 = r[1].z;
	double x3 = r[2].x, y3 = r[2].y, z3 = r[2].z;
	double x4 = r[3].x, y4 = r[3].y, z4 = r[3].z;
	double x5 = r[4].x, y5 = r[4].y, z5 = r[4].z;
	double x6 = r[5].x, y6 = r[5].y, z6 = r[5].z;
	double x7 = r[6].x, y7 = r[6].y, z7 = r[6].z;
	double x8 = r[7].x, y8 = r[7].y, z8 = r[7].z;

	// set up the B-matrix
	double B1[8];
//	double B2[8];
//	double B3[8];

	const double f12 = 1.0/12.0;

	B1[0] = f12*(y2*((z6-z3)-(z4-z5))+y3*(z2-z4)+y4*((z3-z8)-(z5-z2))+y5*((z8-z6)-(z2-z4))+y6*(z5-z2)+y8*(z4-z5));
//	B2[0] = f12*(z2*((x6-x3)-(x4-x5))+z3*(x2-x4)+z4*((x3-x8)-(x5-x2))+z5*((x8-x6)-(x2-x4))+z6*(x5-x2)+z8*(x4-x5));
//	B3[0] = f12*(x2*((y6-y3)-(y4-y5))+x3*(y2-y4)+x4*((y3-y8)-(y5-y2))+x5*((y8-y6)-(y2-y4))+x6*(y5-y2)+x8*(y4-y5));

	B1[1] = f12*(y3*((z7-z4)-(z1-z6))+y4*(z3-z1)+y1*((z4-z5)-(z6-z3))+y6*((z5-z7)-(z3-z1))+y7*(z6-z3)+y5*(z1-z6));
//	B2[1] = f12*(z3*((x7-x4)-(x1-x6))+z4*(x3-x1)+z1*((x4-x5)-(x6-x3))+z6*((x5-x7)-(x3-x1))+z7*(x6-x3)+z5*(x1-x6));
//	B3[1] = f12*(x3*((y7-y4)-(y1-y6))+x4*(y3-y1)+x1*((y4-y5)-(y6-y3))+x6*((y5-y7)-(y3-y1))+x7*(y6-y3)+x5*(y1-y6));

	B1[2] = f12*(y4*((z8-z1)-(z2-z7))+y1*(z4-z2)+y2*((z1-z6)-(z7-z4))+y7*((z6-z8)-(z4-z2))+y8*(z7-z4)+y6*(z2-z7));
//	B2[2] = f12*(z4*((x8-x1)-(x2-x7))+z1*(x4-x2)+z2*((x1-x6)-(x7-x4))+z7*((x6-x8)-(x4-x2))+z8*(x7-x4)+z6*(x2-x7));
//	B3[2] = f12*(x4*((y8-y1)-(y2-y7))+x1*(y4-y2)+x2*((y1-y6)-(y7-y4))+x7*((y6-y8)-(y4-y2))+x8*(y7-y4)+x6*(y2-y7));

	B1[3] = f12*(y1*((z5-z2)-(z3-z8))+y2*(z1-z3)+y3*((z2-z7)-(z8-z1))+y8*((z7-z5)-(z1-z3))+y5*(z8-z1)+y7*(z3-z8));
//	B2[3] = f12*(z1*((x5-x2)-(x3-x8))+z2*(x1-x3)+z3*((x2-x7)-(x8-x1))+z8*((x7-x5)-(x1-x3))+z5*(x8-x1)+z7*(x3-x8));
//	B3[3] = f12*(x1*((y5-y2)-(y3-y8))+x2*(y1-y3)+x3*((y2-y7)-(y8-y1))+x8*((y7-y5)-(y1-y3))+x5*(y8-y1)+x7*(y3-y8));

	B1[4] = f12*(y8*((z4-z7)-(z6-z1))+y7*(z8-z6)+y6*((z7-z2)-(z1-z8))+y1*((z2-z4)-(z8-z6))+y4*(z1-z8)+y2*(z6-z1));
//	B2[4] = f12*(z8*((x4-x7)-(x6-x1))+z7*(x8-x6)+z6*((x7-x2)-(x1-x8))+z1*((x2-x4)-(x8-x6))+z4*(x1-x8)+z2*(x6-x1));
//	B3[4] = f12*(x8*((y4-y7)-(y6-y1))+x7*(y8-y6)+x6*((y7-y2)-(y1-y8))+x1*((y2-y4)-(y8-y6))+x4*(y1-y8)+x2*(y6-y1));

	B1[5] = f12*(y5*((z1-z8)-(z7-z2))+y8*(z5-z7)+y7*((z8-z3)-(z2-z5))+y2*((z3-z1)-(z5-z7))+y1*(z2-z5)+y3*(z7-z2));
//	B2[5] = f12*(z5*((x1-x8)-(x7-x2))+z8*(x5-x7)+z7*((x8-x3)-(x2-x5))+z2*((x3-x1)-(x5-x7))+z1*(x2-x5)+z3*(x7-x2));
//	B3[5] = f12*(x5*((y1-y8)-(y7-y2))+x8*(y5-y7)+x7*((y8-y3)-(y2-y5))+x2*((y3-y1)-(y5-y7))+x1*(y2-y5)+x3*(y7-y2));

	B1[6] = f12*(y6*((z2-z5)-(z8-z3))+y5*(z6-z8)+y8*((z5-z4)-(z3-z6))+y3*((z4-z2)-(z6-z8))+y2*(z3-z6)+y4*(z8-z3));
//	B2[6] = f12*(z6*((x2-x5)-(x8-x3))+z5*(x6-x8)+z8*((x5-x4)-(x3-x6))+z3*((x4-x2)-(x6-x8))+z2*(x3-x6)+z4*(x8-x3));
//	B3[6] = f12*(x6*((y2-y5)-(y8-y3))+x5*(y6-y8)+x8*((y5-y4)-(y3-y6))+x3*((y4-y2)-(y6-y8))+x2*(y3-y6)+x4*(y8-y3));

	B1[7] = f12*(y7*((z3-z6)-(z5-z4))+y6*(z7-z5)+y5*((z6-z1)-(z4-z7))+y4*((z1-z3)-(z7-z5))+y3*(z4-z7)+y1*(z5-z4));
//	B2[7] = f12*(z7*((x3-x6)-(x5-x4))+z6*(x7-x5)+z5*((x6-x1)-(x4-x7))+z4*((x1-x3)-(x7-x5))+z3*(x4-x7)+z1*(x5-x4));
//	B3[7] = f12*(x7*((y3-y6)-(y5-y4))+x6*(y7-y5)+x5*((y6-y1)-(y4-y7))+x4*((y1-y3)-(y7-y5))+x3*(y4-y7)+x1*(y5-y4));

	// calculate the volume V= xi*B1[i] = yi*B2[i] = zi*B3[i] (sum over i)
	return (x1*B1[0]+x2*B1[1]+x3*B1[2]+x4*B1[3]+x5*B1[4]+x6*B1[5]+x7*B1[6]+x8*B1[7]);
}

//-----------------------------------------------------------------------------
//! Assigns a material ID to the entire mesh

void FEMesh::SetMatID(int n)
{
	int i;
	int N = SolidElements();
	for (i=0; i<N; ++i) SolidElement(i).SetMatID(n);
	N = ShellElements();
	for (i=0; i<N; ++i) ShellElement(i).SetMatID(n);
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

		// write mesh item counts
		int nn   = Nodes();
		int nbel = SolidElements();
		int nsel = ShellElements();
		ar << nn << nbel << nsel;

		// write nodal data
		for (i=0; i<nn; ++i) ar.write(&Node(i), sizeof(FENode), 1);

		// write solid element data
		int nmat;
		for (i=0; i<nbel; ++i)
		{
			FESolidElement& el = SolidElement(i);
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
			FEShellElement& el = ShellElement(i);
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
		int nn, nbel, nsel;
		ar >> nn >> nbel >> nsel;

		// allocate storage for mesh data
		Create(nn, nbel, nsel);

		// read nodal data
		for (i=0; i<nn; ++i) ar.read(&Node(i), sizeof(FENode), 1);

		// read solid element data
		for (i=0; i<nbel; ++i)
		{
			FESolidElement& el = SolidElement(i);
			ar >> n;

			switch (n)
			{
			case FE_HEX:
			case FE_RIHEX:
			case FE_PENTA:
			case FE_TET:
				{
					el.SetType(n);

					ar >> mat; el.SetMatID(mat);
					ar >> el.m_nrigid;
					ar >> el.m_nID;
					ar >> el.m_node;

					ar >> el.m_eJ;
					ar >> el.m_ep;
					ar >> el.m_Lk;
				}
				break;
			}
		}

		// read shell element data
		for (i=0; i<nsel; ++i)
		{
			FEShellElement& el = ShellElement(i);
			ar >> n;

			switch (n)
			{
			case FE_SHELL_QUAD:
			case FE_SHELL_TRI:
				{
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
				break;
			}
		}

		// read bounding box data
		FE_BOUNDING_BOX& box = GetBoundingBox();
		ar >> box.r0 >> box.r1;
	}
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
	int i;
	// search the solid elements
	for (i=0; i<m_Elem.size(); ++i)
		if (m_Elem[i].m_nID == nid) return &m_Elem[i];

	// now do the shells
	for (i=0; i<m_Shell.size(); ++i)
		if (m_Shell[i].m_nID == nid) return &m_Shell[i];

	// we could not find it
	return 0;
}
