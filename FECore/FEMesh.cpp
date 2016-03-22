// FEMesh.cpp: implementation of the FEMesh class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMesh.h"
#include "FEException.h"
#include "FEDiscreteDomain.h"
#include "FETrussDomain.h"
#include "FEShellDomain.h"
#include "FEFergusonShellDomain.h"
#include "FESolidDomain.h"
#include "FEDomain2D.h"
#include "FEMaterial.h"
#include "FEModel.h"
#include "log.h"
#include "DOFS.h"
#include "FEElemElemList.h"
#include "FEElementList.h"

//=============================================================================
// FENode
//-----------------------------------------------------------------------------
FENode::FENode()
{
	// exclude flag (true if the node should not be part of the analysis.
	// For instance, if it is isolated).
	m_bexclude = false;

	// shell flag
	m_bshell = false;

	// rigid body data
	m_rid = -1;

	// default ID
	m_nID = -1;
}

//-----------------------------------------------------------------------------
void FENode::SetDOFS(int n)
{
	// initialize dof stuff
    m_ID.assign(n, DOF_FIXED);
	m_BC.assign(n, DOF_OPEN );
	m_val.assign(n, 0.0);
}

//-----------------------------------------------------------------------------
FENode::FENode(const FENode& n)
{
	m_r0 = n.m_r0;
	m_rt = n.m_rt;
	m_at = n.m_at;
	m_rp = n.m_rp;
	m_vp = n.m_vp;
	m_ap = n.m_ap;
	m_Fr = n.m_Fr;

	m_rid = n.m_rid;
	m_bshell = n.m_bshell;
	m_bexclude = n.m_bexclude;

	m_ID = n.m_ID;
	m_BC = n.m_BC;
	m_val = n.m_val;
}

//-----------------------------------------------------------------------------
FENode& FENode::operator = (const FENode& n)
{
	m_r0 = n.m_r0;
	m_rt = n.m_rt;
	m_at = n.m_at;
	m_rp = n.m_rp;
	m_vp = n.m_vp;
	m_ap = n.m_ap;
	m_Fr = n.m_Fr;

	m_rid = n.m_rid;
	m_bshell = n.m_bshell;
	m_bexclude = n.m_bexclude;

	m_ID = n.m_ID;
	m_BC = n.m_BC;
	m_val = n.m_val;

	return (*this);
}

//=============================================================================
// FENodeSet
//-----------------------------------------------------------------------------
FENodeSet::FENodeSet(FEMesh* pm) : m_pmesh(pm), m_nID(-1)
{ 
	m_szname[0] = 0; 
}

//-----------------------------------------------------------------------------
void FENodeSet::create(int n)
{
	assert(n);
	m_Node.resize(n);
}

//-----------------------------------------------------------------------------
void FENodeSet::add(int id)
{
	m_Node.push_back(id);
}

//-----------------------------------------------------------------------------
void FENodeSet::add(const FENodeSet& ns)
{
	int n0 = (int) m_Node.size();
	int n1 = ns.size();
	int N = n0 + n1;
	m_Node.resize(N);
	for (int i=0; i<n1; ++i) m_Node[n0 + i] = ns[i];
}

//-----------------------------------------------------------------------------
void FENodeSet::SetName(const char* sz)
{
	strcpy(m_szname, sz); 
}

//=============================================================================
// FEDiscreteSet
//-----------------------------------------------------------------------------

FEDiscreteSet::FEDiscreteSet(FEMesh* pm) : m_pmesh(pm)
{

}

//-----------------------------------------------------------------------------
void FEDiscreteSet::create(int n)
{
	m_pair.resize(n);
}

//-----------------------------------------------------------------------------
void FEDiscreteSet::add(int n0, int n1)
{
	NodePair p = {n0, n1};
	m_pair.push_back(p);
}

//-----------------------------------------------------------------------------
void FEDiscreteSet::SetName(const char* sz)
{
	strcpy(m_szname, sz); 
}

//=============================================================================
// FEFacetSet
//-----------------------------------------------------------------------------
FEFacetSet::FEFacetSet()
{
	m_szname[0] = 0;
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
void FEFacetSet::SetName(const char* sz)
{
	strcpy(m_szname, sz); 
}

//=============================================================================
// FESegmentSet
//-----------------------------------------------------------------------------
FESegmentSet::FESegmentSet()
{
	m_szname[0] = 0;
}

//-----------------------------------------------------------------------------
void FESegmentSet::Create(int n)
{
	m_Seg.resize(n);
}

//-----------------------------------------------------------------------------
FESegmentSet::SEGMENT& FESegmentSet::Segment(int i)
{
	return m_Seg[i];
}

//-----------------------------------------------------------------------------
void FESegmentSet::SetName(const char* sz)
{
	strcpy(m_szname, sz); 
}

//=============================================================================
// FEElementSet
//-----------------------------------------------------------------------------
FEElementSet::FEElementSet(FEMesh* pm) : m_pmesh(pm)
{
	m_szname[0] = 0;
}

//-----------------------------------------------------------------------------
void FEElementSet::create(int n)
{
	assert(n);
	m_Elem.resize(n);
}

//-----------------------------------------------------------------------------
void FEElementSet::SetName(const char* sz)
{
	strcpy(m_szname, sz); 
}

//=============================================================================
// FEMesh
//-----------------------------------------------------------------------------
FEMesh::FEMesh()
{
}

//-----------------------------------------------------------------------------
FEMesh::~FEMesh()
{
	Clear();
}

//-----------------------------------------------------------------------------
void FEMesh::Serialize(DumpStream& ar)
{
	if (ar.IsShallow())
	{
 		// stream nodal data
		if (ar.IsSaving())
		{
			int NN = (int) m_Node.size();
			for (int i=0; i<NN; ++i)
			{
				FENode& nd = m_Node[i];
				ar << nd.m_r0;
				ar << nd.m_rt << nd.m_at;
				ar << nd.m_rp << nd.m_vp << nd.m_ap;
				ar << nd.m_Fr;
				ar << nd.m_val;
			}
		}
		else
		{
			int NN = (int) m_Node.size();
			for (int i=0; i<NN; ++i)
			{
				FENode& nd = m_Node[i];
				ar >> nd.m_r0;
				ar >> nd.m_rt >> nd.m_at;
				ar >> nd.m_rp >> nd.m_vp >> nd.m_ap;
				ar >> nd.m_Fr;
				ar >> nd.m_val;
			}
		}

		// stream domain data
		int ND = Domains();
		for (int i=0; i<ND; ++i)
		{
			FEDomain& dom = Domain(i);
			dom.Serialize(ar);
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			// write nodal data
			int nn = Nodes();
			ar << nn;
			for (int i=0; i<nn; ++i)
			{
				FENode& node = Node(i);
				ar << node.m_ap;
				ar << node.m_at;
				ar << node.m_bshell;
				ar << node.m_bexclude;
				ar << node.m_Fr;
				ar << node.m_ID;
				ar << node.m_BC;
				ar << node.m_r0;
				ar << node.m_rid;
				ar << node.m_rp;
				ar << node.m_rt;
				ar << node.m_vp;
				ar << node.m_val;
			}

			// write domain data
			int ND = Domains();
			ar << ND;
			for (int i=0; i<ND; ++i)
			{
				FEDomain& d = Domain(i);
				ar << d.GetMaterial()->GetID();
				ar << d.GetTypeStr() << d.Elements();
				d.Serialize(ar);
			}
		}
		else
		{
			FEModel& fem = ar.GetFEModel();
			FECoreKernel& febio = FECoreKernel::GetInstance();

			// read nodal data
			int nn;
			ar >> nn;
			CreateNodes(nn);
			for (int i=0; i<nn; ++i)
			{
				FENode& node = Node(i);
				ar >> node.m_ap;
				ar >> node.m_at;
				ar >> node.m_bshell;
				ar >> node.m_bexclude;
				ar >> node.m_Fr;
				ar >> node.m_ID;
				ar >> node.m_BC;
				ar >> node.m_r0;
				ar >> node.m_rid;
				ar >> node.m_rp;
				ar >> node.m_rt;
				ar >> node.m_vp;
				ar >> node.m_val;
			}

			// read domain data
			int ND, ne;
			ar >> ND;
			char sz[256] = {0};
			for (int i=0; i<ND; ++i)
			{
				int nmat;
				ar >> nmat;
				FEMaterial* pm = fem.FindMaterial(nmat);
				assert(pm);

				ar >> sz >> ne;
				FEDomain* pd = fecore_new<FEDomain>(FEDOMAIN_ID, sz, &fem);
				assert(pd);
				pd->SetMaterial(pm);
				pd->create(ne);
				pd->Serialize(ar);

				AddDomain(pd);
			}

			UpdateBox();
		}
	}
}

//-----------------------------------------------------------------------------
//  Allocates storage for mesh data.
//
void FEMesh::CreateNodes(int nodes)
{
	assert(nodes);
	m_Node.resize(nodes);

	// set the default node IDs
	for (int i=0; i<nodes; ++i) Node(i).SetID(i+1);
}

//-----------------------------------------------------------------------------
// Make more room for nodes
void FEMesh::AddNodes(int nodes)
{
	assert(nodes);
	int N0 = (int) m_Node.size();

	// get the ID of the last node
	// (It is assumed that nodes are sorted according their ID
	//  so the last node should have the highest ID)
	int n0 = 1;
	if (N0 > 0) n0 = m_Node[N0-1].GetID() + 1;

	m_Node.resize(N0 + nodes);
	for (int i=0; i<nodes; ++i) m_Node[i+N0].SetID(n0+i);
}

//-----------------------------------------------------------------------------
void FEMesh::SetDOFS(int n)
{
	int NN = Nodes();
	for (int i=0; i<NN; ++i) m_Node[i].SetDOFS(n);
}

//-----------------------------------------------------------------------------
//! Return the total number elements
int FEMesh::Elements() const
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i) 
	{
		N += m_Domain[i]->Elements();
	}
	return N;
}

//-----------------------------------------------------------------------------
//! Return the total number of elements of a specific domain type
int FEMesh::Elements(int ndom_type) const
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i) 
	{
		FEDomain& dom = *m_Domain[i];
		if (dom.Class() == ndom_type) N += m_Domain[i]->Elements();
	}
	return N;
}

//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
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

	// See if there are any isolated nodes
	// Exclude them from the analysis
	int ni = 0;
	for (i=0; i<N; ++i)
		if (val[i] == 0)
		{
			++ni;
			FENode& node = Node(i);
			node.m_bexclude = true;
		}

	return ni;
}

//-----------------------------------------------------------------------------
//! Calculate all shell normals (i.e. the shell directors).
//! And find shell nodes
void FEMesh::InitShells()
{
	// zero initial directors for shell nodes
	int NN = Nodes();
	vector<vec3d> D(NN, vec3d(0,0,0));

	// loop over all domains
	for (int nd = 0; nd < Domains(); ++nd)
	{
		// Calculate the shell directors as the local node normals
		if (Domain(nd).Class() == FE_DOMAIN_SHELL)
		{
			FEShellDomain& sd = static_cast<FEShellDomain&>(Domain(nd));
			vec3d r0[FEElement::MAX_NODES];
			for (int i=0; i<sd.Elements(); ++i)
			{
				FEShellElement& el = sd.Element(i);

				int n = el.Nodes();
				int* en = &el.m_node[0];

				// get the nodes
				for (int j=0; j<n; ++j) r0[j] = Node(en[j]).m_r0;

				for (int j=0; j<n; ++j)
				{
					int m0 = j;
					int m1 = (j+1)%n;
					int m2 = (j==0? n-1: j-1);

					vec3d a = r0[m0];
					vec3d b = r0[m1];
					vec3d c = r0[m2];

					D[en[m0]] += (b-a)^(c-a);
				}
			}
		}
        else if (Domain(nd).Class() == FE_DOMAIN_FERGUSON)
        {
            FEFergusonShellDomain& sd = static_cast<FEFergusonShellDomain&>(Domain(nd));
            vec3d r0[FEElement::MAX_NODES];
            for (int i=0; i<sd.Elements(); ++i)
            {
                FEFergusonShellElement& el = sd.Element(i);
                
                int n = el.Nodes();
                int* en = &el.m_node[0];
                
                // get the nodes
                for (int j=0; j<n; ++j) r0[j] = Node(en[j]).m_r0;
                
                for (int j=0; j<n; ++j)
                {
                    int m0 = j;
                    int m1 = (j+1)%n;
                    int m2 = (j==0? n-1: j-1);
                    
                    vec3d a = r0[m0];
                    vec3d b = r0[m1];
                    vec3d c = r0[m2];
                    
                    D[en[m0]] += (b-a)^(c-a);
                }
            }
        }
	}

	// make sure we start with unit directors
	for (int i=0; i<NN; ++i) D[i].unit();

	// assign directors to shells 
	for (int nd = 0; nd < Domains(); ++nd)
	{
		// Calculate the shell directors as the local node normals
		if (Domain(nd).Class() == FE_DOMAIN_SHELL)
		{
			FEShellDomain& sd = static_cast<FEShellDomain&>(Domain(nd));
			for (int i=0; i<sd.Elements(); ++i)
			{
				FEShellElement& el = sd.Element(i);
				int ne = el.Nodes();
				for (int j=0; j<ne; ++j) el.m_D0[j] = D[el.m_node[j]]*el.m_h0[j];
			}
		}
        else if (Domain(nd).Class() == FE_DOMAIN_FERGUSON)
        {
            FEFergusonShellDomain& sd = static_cast<FEFergusonShellDomain&>(Domain(nd));
            for (int i=0; i<sd.Elements(); ++i)
            {
                FEFergusonShellElement& el = sd.Element(i);
                int ne = el.Nodes();
                for (int j=0; j<ne; ++j) el.m_D0[j] = D[el.m_node[j]]*el.m_h0[j];
            }
        }
	}

	// Find the nodes that are on a non-rigid shell. 
	// These nodes will be assigned rotational degrees of freedom
	// TODO: Perhaps I should let the domains do this instead
	for (int i=0; i<Nodes(); ++i) Node(i).m_bshell = false;
	for (int nd = 0; nd<Domains(); ++nd)
	{
		FEDomain& dom = Domain(nd);
		if ((dom.Class() == FE_DOMAIN_SHELL) || (dom.Class() == FE_DOMAIN_FERGUSON))
		{
			FEMaterial* pmat = dom.GetMaterial();
			if (pmat->IsRigid() == false)
			{
				int N = dom.Elements();
				for (int i=0; i<N; ++i)
				{
					FEElement& el = dom.ElementRef(i);
					int n = el.Nodes();
					for (int j=0; j<n; ++j) Node(el.m_node[j]).m_bshell = true;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Does one-time initialization of the Mesh data. Call FEMesh::Reset for resetting 
//! the mesh data.
bool FEMesh::Init()
{
	// find and remove isolated vertices
	int ni = RemoveIsolatedVertices();
	if (ni != 0) 
	{
		if (ni == 1)
			felog.printbox("WARNING", "%d isolated vertex removed.", ni);
		else
			felog.printbox("WARNING", "%d isolated vertices removed.", ni);
	}

	// Initialize shell normals (i.e. directors)
	// NOTE: we do this before we check for inverted elements since the jacobian of a shell
	//       depends on its normal.
	InitShells();

	// reset data
	// TODO: Not sure why this is here
	Reset();

	// All done
	return true;
}

//-----------------------------------------------------------------------------
void FEMesh::Clear()
{
	m_Node.clear();
	for (size_t i=0; i<m_Domain.size (); ++i) delete m_Domain [i];

	// TODO: Surfaces are currently managed by the classes that use them so don't delete them
//	for (size_t i=0; i<m_Surf.size   (); ++i) delete m_Surf   [i];

	for (size_t i=0; i<m_NodeSet.size(); ++i) delete m_NodeSet[i];
	for (size_t i=0; i<m_LineSet.size(); ++i) delete m_LineSet[i];
	for (size_t i=0; i<m_FaceSet.size(); ++i) delete m_FaceSet[i];
	for (size_t i=0; i<m_ElemSet.size  (); ++i) delete m_ElemSet[i];
	m_Domain.clear();
	m_Surf.clear();
	m_NodeSet.clear();
	m_LineSet.clear();
	m_FaceSet.clear();
	m_ElemSet.clear();
	m_NEL.Clear();
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
		node.m_vp = vec3d(0,0,0);
		node.m_ap = node.m_at = vec3d(0,0,0);

        node.m_Fr = vec3d(0,0,0);

		// reset ID arrays
		int ndof = node.m_ID.size();
		for (int i=0; i<ndof; ++i) 
		{
			node.m_ID[i] = DOF_FIXED;
			node.m_BC[i] = DOF_OPEN;
			node.set(i, 0.0);
		}
	}

	// update the mesh
	UpdateBox();

	// reset domain data
	for (int n=0; n<(int) m_Domain.size(); ++n) m_Domain[n]->Reset();
}

//-----------------------------------------------------------------------------
//! This function calculates the (initial) volume of an element. In some case, the volume
//! may only be approximate.
double FEMesh::ElementVolume(FEElement &el)
{
	double V = 0;
	switch (el.Class())
	{
	case FE_ELEM_SOLID: V = SolidElementVolume(static_cast<FESolidElement&>(el)); break;
	case FE_ELEM_SHELL: V = ShellElementVolume(static_cast<FEShellElement&>(el)); break;
	}
	return V;
}

//-----------------------------------------------------------------------------
//! \todo Replace this with what FEBio 1.x does.
double FEMesh::SolidElementVolume(FESolidElement& el)
{
	int i;
	vec3d r0[FEElement::MAX_NODES];

	int neln = el.Nodes();
	for (i=0; i<neln; ++i) r0[i] = Node(el.m_node[i]).m_r0;

	int nint = el.GaussPoints();
	double *w = el.GaussWeights();
	double V = 0;
	for (int n=0; n<nint; ++n) 
	{
		// shape function derivatives
		double* Grn = el.Gr(n);
		double* Gsn = el.Gs(n);
		double* Gtn = el.Gt(n);

		// jacobian matrix
		double J[3][3] = {0};
		for (i=0; i<neln; ++i)
		{
			const double& Gri = Grn[i];
			const double& Gsi = Gsn[i];
			const double& Gti = Gtn[i];
			
			const double& x = r0[i].x;
			const double& y = r0[i].y;
			const double& z = r0[i].z;
			
			J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
			J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
			J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
		}
			
		// calculate the determinant
		double detJ0 =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
					+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
					+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

		V += detJ0*w[n];
	}

	return V;
}

//-----------------------------------------------------------------------------
//! \todo Replace this with what FEBio 1.x does.
double FEMesh::ShellElementVolume(FEShellElement& el)
{
	int i;
	int neln = el.Nodes();

	// initial nodal coordinates and directors
	vec3d r0[FEElement::MAX_NODES], D0[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		r0[i] = Node(el.m_node[i]).m_r0;
		D0[i] = el.m_D0[i];
	}

	int nint = el.GaussPoints();
	double *w = el.GaussWeights();
	double V = 0;
    vec3d g[3];
	for (int n=0; n<nint; ++n)
	{
		// jacobian matrix
		double eta = el.gt(n);
        
        double* Mr = el.Hr(n);
        double* Ms = el.Hs(n);
        double* M  = el.H(n);
        
        // evaluate covariant basis vectors
        g[0] = g[1] = g[2] = vec3d(0,0,0);
        for (i=0; i<neln; ++i)
        {
            g[0] += (r0[i] + D0[i]*eta/2)*Mr[i];
            g[1] += (r0[i] + D0[i]*eta/2)*Ms[i];
            g[2] += D0[i]*(M[i]/2);
        }
				
        mat3d J = mat3d(g[0].x, g[1].x, g[2].x,
                        g[0].y, g[1].y, g[2].y,
                        g[0].z, g[1].z, g[2].z);
        
		// calculate the determinant
        double detJ0 = J.det();

		V += detJ0*w[n];
	}

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
	for (size_t i=0; i<m_NodeSet.size(); ++i) if (strcmp(m_NodeSet[i]->GetName(), szname) == 0) return m_NodeSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a segment set set by name

FESegmentSet* FEMesh::FindSegmentSet(const char* szname)
{
	for (size_t i=0; i<m_LineSet.size(); ++i) if (strcmp(m_LineSet[i]->GetName(), szname) == 0) return m_LineSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a discrete element set set by name

FEDiscreteSet* FEMesh::FindDiscreteSet(const char* szname)
{
	for (size_t i=0; i<m_DiscSet.size(); ++i) if (strcmp(m_DiscSet[i]->GetName(), szname) == 0) return m_DiscSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a element set by name

FEElementSet* FEMesh::FindElementSet(const char* szname)
{
	for (size_t i=0; i<m_ElemSet.size(); ++i) if (strcmp(m_ElemSet[i]->GetName(), szname) == 0) return m_ElemSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
int FEMesh::Faces(FEElement& el)
{
	switch (el.Type())
	{
	case FE_HEX8G8:
	case FE_HEX8RI:
	case FE_HEX8G1:
	case FE_HEX20G27: 
	case FE_HEX27G27: return 6;
	case FE_PENTA6G6: return 5;
	case FE_TET4G4:
	case FE_TET10G4:
	case FE_TET10G8:
	case FE_TET10GL11:
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
	case FE_TET4G1: return 4;
	case FE_SHELL_QUAD:
    case FE_SHELL_QUAD8:
	case FE_SHELL_TRI:
    case FE_SHELL_TRI6: return 1;
	default:
		assert(false);
	}

	return 0;
}

//-----------------------------------------------------------------------------
//! This function returns the face connectivity from a certain element

int FEMesh::GetFace(FEElement& el, int n, int* nf)
{
	int nn = -1;
	int* en = &el.m_node[0];
	switch (el.Type())
	{
	case FE_HEX8G8:
	case FE_HEX8RI:
	case FE_HEX8G1:
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
	case FE_PENTA6G6:
		switch(n)
		{
		case 0: nn = 4; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; nf[3] = en[3]; break;
		case 1: nn = 4; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[5]; nf[3] = en[4]; break;
		case 2: nn = 4; nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[5]; nf[3] = en[2]; break;
		case 3: nn = 3; nf[0] = en[0]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[1]; break;
		case 4: nn = 3; nf[0] = en[3]; nf[1] = en[4]; nf[2] = en[5]; nf[3] = en[5]; break;
		}
		break;
	case FE_TET4G4:
	case FE_TET4G1:
		nn = 3;
		switch (n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = nf[3] = en[3]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = nf[3] = en[3]; break;
		case 2: nf[0] = en[0]; nf[1] = en[3]; nf[2] = nf[3] = en[2]; break;
		case 3: nf[0] = en[0]; nf[1] = en[2]; nf[2] = nf[3] = en[1]; break;
		}
		break;
	case FE_TET10G4:
	case FE_TET10G8:
	case FE_TET10GL11:
		nn = 6;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[7]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[8]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[9]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[5]; nf[4] = en[4]; nf[5] = en[6]; break;
		}
		break;
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
		nn = 7;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[7]; nf[6] = en[11]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[8]; nf[6] = en[12]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[9]; nf[6] = en[13]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[5]; nf[4] = en[4]; nf[5] = en[6]; nf[6] = en[10]; break;
		}
		break;
	case FE_HEX20G27:
		nn = 8;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[ 8]; nf[5] = en[17]; nf[6] = en[12]; nf[7] = en[16]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; nf[4] = en[ 9]; nf[5] = en[18]; nf[6] = en[13]; nf[7] = en[17]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; nf[4] = en[10]; nf[5] = en[19]; nf[6] = en[14]; nf[7] = en[18]; break;
		case 3: nf[0] = en[0]; nf[1] = en[4]; nf[2] = en[7]; nf[3] = en[3]; nf[4] = en[16]; nf[5] = en[15]; nf[6] = en[19]; nf[7] = en[11]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; nf[4] = en[11]; nf[5] = en[10]; nf[6] = en[ 9]; nf[7] = en[ 8]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; nf[4] = en[12]; nf[5] = en[13]; nf[6] = en[14]; nf[7] = en[15]; break;
		}
		break;
	case FE_HEX27G27:
		nn = 9;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[ 8]; nf[5] = en[17]; nf[6] = en[12]; nf[7] = en[16]; nf[8] = en[20]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; nf[4] = en[ 9]; nf[5] = en[18]; nf[6] = en[13]; nf[7] = en[17]; nf[8] = en[21]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; nf[4] = en[10]; nf[5] = en[19]; nf[6] = en[14]; nf[7] = en[18]; nf[8] = en[22]; break;
		case 3: nf[0] = en[0]; nf[1] = en[4]; nf[2] = en[7]; nf[3] = en[3]; nf[4] = en[16]; nf[5] = en[15]; nf[6] = en[19]; nf[7] = en[11]; nf[8] = en[23]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; nf[4] = en[11]; nf[5] = en[10]; nf[6] = en[ 9]; nf[7] = en[ 8]; nf[8] = en[24]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; nf[4] = en[12]; nf[5] = en[13]; nf[6] = en[14]; nf[7] = en[15]; nf[8] = en[25]; break;
		}
		break;
	case FE_SHELL_QUAD:
		nn = 4;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3];
		break;
    case FE_SHELL_QUAD8:
        nn = 8;
        nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6]; nf[7] = en[7];
        break;
	case FE_SHELL_TRI:
		nn = 3;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2];
		break;
    case FE_SHELL_TRI6:
        nn = 6;
        nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5];
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
		if (m_Domain[i]->Class() == FE_DOMAIN_SOLID)
		{
			FESolidDomain& bd = static_cast<FESolidDomain&>(*m_Domain[i]);
			FESolidElement* pe = bd.FindElement(y, r);
			if (pe) return pe;
		}
	}
	return 0;
}


//-----------------------------------------------------------------------------
//! This function finds all the domains that have a certain material
void FEMesh::DomainListFromMaterial(vector<int>& lmat, vector<int>& ldom)
{
	// make sure the list is empty
	if (ldom.empty() == false) ldom.clear();

	// loop over all domains
	int ND = (int) m_Domain.size();
	int NM = (int) lmat.size();
	for (int i=0; i<ND; ++i)
	{
		FEDomain& di = *m_Domain[i];
		int dmat = di.GetMaterial()->GetID();
		for (int j=0; j<NM; ++j)
		{
			if (dmat == lmat[j])
			{
				ldom.push_back(i);
				break;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the surface representing the element boundaries
//! boutside : include all exterior facets
//! binside  : include all interior facets
FESurface* FEMesh::ElementBoundarySurface(bool boutside, bool binside)
{
	if ((boutside == false) && (binside == false)) return 0;

	// create the element neighbor list
	FEElemElemList EEL;
	EEL.Create(this);

	// get the number of elements in this mesh
	int NE = Elements();

	// count the number of facets we have to create
	int NF = 0;
	FEElementList EL(*this);
	FEElementList::iterator it = EL.begin();
	for (int i=0; i<NE; ++i, ++it)
	{
		FEElement& el = *it;
		int nf = Faces(el);
		for (int j=0; j<nf; ++j)
		{
			FEElement* pen = EEL.Neighbor(i, j);
			if ((pen == 0) && boutside) ++NF;
			if ((pen != 0) && binside ) ++NF;
		}
	}
	// create the surface
	FESurface* ps = new FESurface(this);
	if (NF == 0) return 0;
	ps->create(NF);

	// build the surface elements
	int face[FEElement::MAX_NODES];
	NF = 0;
	it = EL.begin();
	for (int i=0; i<NE; ++i, ++it)
	{
		FEElement& el = *it;
		int nf = Faces(el);
		for (int j=0; j<nf; ++j)
		{
			FEElement* pen = EEL.Neighbor(i, j);
			if (((pen == 0) && boutside)||
				((pen != 0) && binside ))
			{
				FESurfaceElement& se = ps->Element(NF++);
				GetFace(el, j, face);

				switch (el.Type())
				{
				case FE_HEX8G8:
				case FE_SHELL_QUAD:
					se.SetType(FE_QUAD4G4); 
					break;
                case FE_SHELL_QUAD8:
                    se.SetType(FE_QUAD8G9);
                    break;
				case FE_TET4G1: 
				case FE_SHELL_TRI:
					se.SetType(FE_TRI3G1); 
					break;
                case FE_SHELL_TRI6:
                    se.SetType(FE_TRI6G7);
                    break;
				}
				
				se.m_nelem = el.GetID();
				
				int nn = se.Nodes();
				for (int k=0; k<nn; ++k)
				{
					se.m_node[k] = face[k];
				}
			}
		}
	}

	// initialize the surface. 
	// This will set the local surface element ID's and also set the m_nelem IDs.
	ps->Init();

	// all done
	return ps;
}
