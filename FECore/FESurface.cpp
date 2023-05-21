/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FESurface.h"
#include "FEMesh.h"
#include "FESolidDomain.h"
#include "FEElemElemList.h"
#include "DumpStream.h"
#include "matrix.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
FESurface::FESurface(FEModel* fem) : FEMeshPartition(FE_DOMAIN_SURFACE, fem)
{
	m_surf = 0;
	m_bitfc = false;
	m_alpha = 1;
	m_bshellb = false;
}

//-----------------------------------------------------------------------------
FESurface::~FESurface()
{

}

//-----------------------------------------------------------------------------
void FESurface::Create(int nsize, int elemType)
{
	m_el.resize(nsize);
	for (int i = 0; i < nsize; ++i)
	{
		FESurfaceElement& el = m_el[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
		el.m_elem[0] = nullptr;
		el.m_elem[1] = nullptr;
	}

	if (elemType != -1)
	{
		for (int i = 0; i < nsize; ++i) m_el[i].SetType(elemType);
		CreateMaterialPointData();
	}
}

//-----------------------------------------------------------------------------
void FESurface::Create(const FEFacetSet& set)
{
	if (m_surf == 0) m_surf = const_cast<FEFacetSet*>(&set);
	assert(m_surf == &set);

	FEMesh& m = *GetMesh();

	// count nr of faces
	int faces = set.Faces();

	// allocate storage for faces
	Create(faces);

	// read faces
	for (int i = 0; i<faces; ++i)
	{
		FESurfaceElement& el = Element(i);
		const FEFacetSet::FACET& fi = set.Face(i);

		if (fi.ntype == 4) el.SetType(FE_QUAD4G4);
		else if (fi.ntype == 3) el.SetType(FE_TRI3G1);
		else if (fi.ntype == 6) el.SetType(FE_TRI6G3);
		else if (fi.ntype == 7) el.SetType(FE_TRI7G4);
		else if (fi.ntype == 8) el.SetType(FE_QUAD8G9);
		else if (fi.ntype == 9) el.SetType(FE_QUAD9G9);
		else if (fi.ntype == 10) el.SetType(FE_TRI10G7);
		else assert(false);

		int N = el.Nodes(); assert(N == fi.ntype);
		for (int j = 0; j<N; ++j) el.m_node[j] = fi.node[j];
	}

	// copy the name
	SetName(set.GetName());

	// allocate surface material points
	CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
void FESurface::CreateMaterialPointData()
{
	for (int i = 0; i < Elements(); ++i)
	{
		FESurfaceElement& el = m_el[i];
		int nint = el.GaussPoints();
		el.ClearData();
		for (int n = 0; n < nint; ++n)
		{
			FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(CreateMaterialPoint());
			assert(pt);
			el.SetMaterialPointData(pt, n);
		}
	}
}

//-----------------------------------------------------------------------------
// extract the nodes from this surface
FENodeList FESurface::GetNodeList()
{
	FEMesh* pm = GetMesh();
	FENodeList nset(pm);

	vector<int> tag(pm->Nodes(), 0);
	for (int i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		int ne = el.Nodes();
		for (int j=0; j<ne; ++j)
		{
			if (tag[el.m_node[j]] == 0)
			{
				nset.Add(el.m_node[j]);
				tag[el.m_node[j]] = 1;
			}
		}
	}
	return nset;
}

//-----------------------------------------------------------------------------
//! Get the list of (local) node indices of the boundary nodes
void FESurface::GetBoundaryFlags(std::vector<bool>& boundary) const
{
	FEElemElemList EEL;
	EEL.Create(this);

	boundary.assign(Nodes(), false);
	for (int i = 0; i < Elements(); ++i) {
		const FESurfaceElement& el = Element(i);
		for (int j = 0; j < el.facet_edges(); ++j) {
			FEElement* nel = EEL.Neighbor(i, j);
			if (nel == nullptr) {
				int en[3] = { -1,-1,-1 };
				el.facet_edge(j, en);
				boundary[en[0]] = true;
				boundary[en[1]] = true;
				if (en[2] > -1) boundary[en[2]] = true;
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Create material point data for this surface
FEMaterialPoint* FESurface::CreateMaterialPoint()
{
	return new FESurfaceMaterialPoint;
}

//-----------------------------------------------------------------------------
// update surface data
void FESurface::Update(const FETimeInfo& tp)
{
	ForEachSurfaceElement([=](FESurfaceElement& el) {
		int nint = el.GaussPoints();
		int neln = el.Nodes();

		vec3d rt[FEElement::MAX_NODES];
		NodalCoordinates(el, rt);

		for (int n = 0; n < nint; ++n)
		{
			FESurfaceMaterialPoint& mp = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));

			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);

			mp.dxr = vec3d(0, 0, 0);
			mp.dxs = vec3d(0, 0, 0);
			for (int i = 0; i < neln; ++i)
			{
				mp.dxr += rt[i] * Gr[i];
				mp.dxs += rt[i] * Gs[i];
			}
		}
	});
}

//-----------------------------------------------------------------------------
void FESurface::InitSurface()
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *GetMesh();

	// This array is used to keep tags on each node
	vector<int> tag; tag.assign(mesh.Nodes(), -1);

	// let's find all nodes the surface needs
	int nn = 0;
	int ne = Elements();
	for (int i = 0; i<ne; ++i)
	{
		FESurfaceElement& el = Element(i);
		el.m_lid = i;

		for (int j = 0; j<el.Nodes(); ++j)
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
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		if (tag[i] >= 0)
		{
			m_Node[tag[i]] = i;
		}
	}
}

//-----------------------------------------------------------------------------
//! Initialize surface node data structure
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FESurface::Init()
{
	// make sure that there is a surface defined
	if (Elements() == 0) return false;

	// initialize the surface data
	InitSurface();

	// NOTE: Make sure the mesh has its node-element list initialized
	// otherwise the omp loop below might crash due to race condition.
	GetMesh()->NodeElementList();

	// see if we can find all elements that the faces belong to
	int invalidFacets = 0;
	int ne = Elements();
#pragma omp parallel for reduction(+:invalidFacets)
	for (int i=0; i<ne; ++i)
	{
		FESurfaceElement& el = Element(i);
        if (m_bitfc && (el.m_elem[0] == nullptr)) FindElements(el);
		else if (el.m_elem[0] == nullptr) el.m_elem[0] = FindElement(el);
        //to make sure
        else if (m_bitfc && (el.m_elem[1] == nullptr)) FindElements(el);
		if (el.m_elem[0] == nullptr) { invalidFacets++; }
	}

	if (invalidFacets > 0)
	{
		std::string surfName = GetName();
		if (surfName.empty()) surfName = "(unknown)";
		feLogWarning("The surface \"%s\" has %d invalid facets. \nThe model may not run correctly.", surfName.c_str(), invalidFacets);
	}

	vec3d re[FEElement::MAX_NODES];
	// initialize material points of surface elements
	for (int i = 0; i < Elements(); ++i)
	{
		FESurfaceElement& el = m_el[i];

		NodalCoordinates(el, re);

		int nint = el.GaussPoints();
		int neln = el.Nodes();
		for (int n = 0; n < nint; ++n)
		{
			FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(el.GetMaterialPoint(n));
			if (pt == nullptr) return false;

			// initialize some material point data
			double* H = el.H(n);
			vec3d rn(0, 0, 0);
			for (int j = 0; j < neln; ++j)
			{
				rn += re[j]*H[j];
			}

			pt->m_r0 = rn;

			// calculate initial surface tangents
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);

			vec3d dxr(0, 0, 0), dxs(0, 0, 0);
			for (int i = 0; i < neln; ++i)
			{
				dxr += re[i] * Gr[i];
				dxs += re[i] * Gs[i];
			}
			pt->dxr = dxr;
			pt->dxs = dxs;

			// initialize the other material point data
			pt->Init();
		}
	}

    // allocate node normals and evaluate them in initial configuration
    m_nn.assign(Nodes(), vec3d(0,0,0));
    UpdateNodeNormals();
    
	return true;
}

//-----------------------------------------------------------------------------
//! Find the element that a face belongs to
// TODO: I should be able to speed this up
FEElement* FESurface::FindElement(FESurfaceElement& el)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *GetMesh();
	FENodeElemList& NEL = mesh.NodeElementList();

	int node = el.m_node[0];
	int nval = NEL.Valence(node);
	FEElement** ppe = NEL.ElementList(node);
	for (int i=0; i<nval; ++i)
	{
		FEElement* pe = ppe[i];
		int nfaces = pe->Faces();
		
		int nf[FEElement::MAX_NODES], nn;
		for (int j=0; j<nfaces; ++j)
		{
			nn = pe->GetFace(j, nf);
			if (nn == el.Nodes())
			{
				switch (nn)
				{
				case  3: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) return pe; break;
				case  4: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2]) && el.HasNode(nf[3])) return pe; break;
				case  6: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) return pe; break;
				case  7: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) return pe; break;
				case  8: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2]) && el.HasNode(nf[3])) return pe; break;
				case  9: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2]) && el.HasNode(nf[3])) return pe; break;
				case 10: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) return pe; break;
				default:
					assert(false);
				}
			}
		}
	}

	return nullptr;
}

void FESurface::ForEachSurfaceElement(std::function<void(FESurfaceElement& el)> f)
{
	for (size_t i = 0; i < m_el.size(); ++i) f(m_el[i]);
}

void FESurface::FindElements(FESurfaceElement& el)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *GetMesh();
	FENodeElemList& NEL = mesh.NodeElementList();

	vector<int>& sf = el.m_node;
	int node = el.m_node[0];
	int nval = NEL.Valence(node);
	FEElement** ppe = NEL.ElementList(node);
	for (int i = 0; i < nval; ++i)
	{
		FEElement& sel = *ppe[i];
            
		// check all faces of this solid element
		int nfaces = sel.Faces();
		for (int j = 0; j<nfaces; ++j) 
		{
			int nf[9];
			vec3d g[3];
			int nn = sel.GetFace(j, nf);
                
			int found = 0;
			if (nn == el.Nodes())
			{
                switch (nn)
                {
                    case 3: found = el.HasNodes(nf,3); break;
                    case 4: found = el.HasNodes(nf,4); break;
                    case 6: found = el.HasNodes(nf,3); break;
                    case 7: found = el.HasNodes(nf,3); break;
                    case 8: found = el.HasNodes(nf,4); break;
                    case 9: found = el.HasNodes(nf,4); break;
                    default:
                        assert(false);
                }
            }
            if (found != 0) {
                if (el.m_elem[0] == nullptr) { el.m_elem[0] = &sel; }
                else if (el.m_elem[0] != &sel) el.m_elem[1] = &sel;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! unpack an LM vector from a dof list
void FESurface::UnpackLM(const FESurfaceElement& el, const FEDofList& dofList, vector<int>& lm)
{
	int dofPerNode = dofList.Size();
	int neln = el.Nodes();
	int ndof = neln*dofPerNode;
	lm.assign(ndof, -1);
	for (int j = 0; j < neln; ++j)
	{
		FENode& node = Node(el.m_lnode[j]);
		for (int k = 0; k < dofPerNode; ++k)
		{
			if (dofList[k] >= 0)
				lm[dofPerNode*j + k] = node.m_ID[dofList[k]];
		}
	}
}

//-----------------------------------------------------------------------------
// project onto a triangular face
vec3d project2tri(vec3d* y, vec3d x, double& r, double& s)
{
	// calculate base vectors 
	vec3d e1 = y[1] - y[0];
	vec3d e2 = y[2] - y[0];

	// calculate plane normal
	vec3d n = e1^e2; n.unit();

	// project x onto the plane
	vec3d q = x - n*((x-y[0])*n);

	// set up metric tensor
	double G[2][2];
	G[0][0] = e1*e1;
	G[0][1] = G[1][0] = e1*e2;
	G[1][1] = e2*e2;

	// invert metric tensor
	double D = G[0][0]*G[1][1] - G[0][1]*G[1][0];
	double Gi[2][2];
	Gi[0][0] = G[1][1]/D;
	Gi[1][1] = G[0][0]/D;
	Gi[0][1] = Gi[1][0] = -G[0][1]/D;

	// calculate dual base vectors
	vec3d E1 = e1*Gi[0][0] + e2*Gi[0][1];
	vec3d E2 = e1*Gi[1][0] + e2*Gi[1][1];

	// now we can calculate r and s
	vec3d t = q - y[0];
	r = t*E1;
	s = t*E2;

	return q;
}

//-----------------------------------------------------------------------------
// project onto a quadrilateral surface.
bool project2quad(vec3d* y, vec3d x, double& r, double& s, vec3d& q)
{
	double R[2], u[2], D;
	double gr[4] = {-1, +1, +1, -1};
	double gs[4] = {-1, -1, +1, +1};
	double H[4], Hr[4], Hs[4], Hrs[4];

	int i, j;
	int NMAX = 50, n=0;

	// evaulate scalar products
	double xy[4] = {x*y[0], x*y[1], x*y[2], x*y[3]};
	double yy[4][4];
	yy[0][0] = y[0]*y[0]; yy[1][1] = y[1]*y[1]; yy[2][2] = y[2]*y[2]; yy[3][3] = y[3]*y[3];
	yy[0][1] = yy[1][0] = y[0]*y[1];
	yy[0][2] = yy[2][0] = y[0]*y[2];
	yy[0][3] = yy[3][0] = y[0]*y[3];
	yy[1][2] = yy[2][1] = y[1]*y[2];
	yy[1][3] = yy[3][1] = y[1]*y[3];
	yy[2][3] = yy[3][2] = y[2]*y[3];

	// loop until converged
	bool bconv = false;
	double normu;
	do
	{
		// evaluate shape functions and shape function derivatives.
		for (i=0; i<4; ++i)
		{
			H[i] = 0.25*(1+gr[i]*r)*(1+gs[i]*s);
	
			Hr[i] = 0.25*gr[i]*( 1 + gs[i]*s );
			Hs[i] = 0.25*gs[i]*( 1 + gr[i]*r );

			Hrs[i] = 0.25*gr[i]*gs[i];
		}

		// set up the system of equations
		R[0] = R[1] = 0;
		double A[2][2] = {0};
		for (i=0; i<4; ++i)
		{
			R[0] -= (xy[i])*Hr[i];
			R[1] -= (xy[i])*Hs[i];

			A[0][1] += (xy[i])*Hrs[i];
			A[1][0] += (xy[i])*Hrs[i];

			for (j=0; j<4; ++j)
			{
				double yij = yy[i][j];
				R[0] -= -H[j]*Hr[i]*(yij);
				R[1] -= -H[j]*Hs[i]*(yij);

				A[0][0] -= (yij)*(Hr[i]*Hr[j]);
				A[1][1] -= (yij)*(Hs[i]*Hs[j]);

				A[0][1] -= (yij)*(Hr[i]*Hs[j]+Hrs[i]*H[j]);
				A[1][0] -= (yij)*(Hs[i]*Hr[j]+Hrs[i]*H[j]);
			}
		}
	
		// determinant of A
		D = A[0][0]*A[1][1] - A[0][1]*A[1][0];

		// solve for u = A^(-1)*R
		u[0] = (A[1][1]*R[0] - A[0][1]*R[1])/D;
		u[1] = (A[0][0]*R[1] - A[1][0]*R[0])/D;

		// calculate displacement norm
		normu = u[0]*u[0]+u[1]*u[1];

		// check for convergence
		bconv = ((normu < 1e-10));
		if (!bconv && (n <= NMAX))
		{
			// Don't update if converged otherwise the point q
			// does not correspond with the current values for (r,s)
			r += u[0];
			s += u[1];
			++n;
		}
		else break;
	}
	while (1);

	// evaluate q
	q = y[0]*H[0] + y[1]*H[1] + y[2]*H[2] + y[3]*H[3];

	return bconv;
}

//-----------------------------------------------------------------------------
// project to general surface element.
bool project2surf(FESurfaceElement& el, vec3d* y, vec3d x, double& r, double& s, vec3d& q)
{
	double Q[2], u[2], D;
	const int NE = FEElement::MAX_NODES;
	double H[NE], Hr[NE], Hs[NE], Hrr[NE], Hss[NE], Hrs[NE];

	int i, j;
	int NMAX = 50, n=0;

	// get number of nodes
	int ne = el.Nodes();

	// evaulate scalar products
	double xy[NE];
	double yy[NE][NE];
	for (i=0; i<ne; ++i)
	{
		xy[i] = x*y[i];
		yy[i][i] = y[i]*y[i];
		for (j=i+1; j<ne; ++j)
		{
			yy[i][j] = yy[j][i] = y[i]*y[j];
		}
	}

	// loop until converged
	bool bconv = false;
	double normu;
	do
	{
		// evaluate shape functions and shape function derivatives.
		el.shape_fnc(H, r, s);
		el.shape_deriv(Hr, Hs, r, s);
		el.shape_deriv2(Hrr, Hrs, Hss, r, s);

		// set up the system of equations
		Q[0] = Q[1] = 0;
		double A[2][2] = {0};
		for (i=0; i<ne; ++i)
		{
			Q[0] -= (xy[i])*Hr[i];
			Q[1] -= (xy[i])*Hs[i];

			A[0][0] += (xy[i])*Hrr[i];
			A[0][1] += (xy[i])*Hrs[i];
			A[1][0] += (xy[i])*Hrs[i];
			A[1][1] += (xy[i])*Hss[i];

			for (j=0; j<ne; ++j)
			{
				double yij = yy[i][j];
				Q[0] -= -H[j]*Hr[i]*(yij);
				Q[1] -= -H[j]*Hs[i]*(yij);

				A[0][0] -= (yij)*(H[i]*Hrr[j] + Hr[i]*Hr[j]);
				A[1][1] -= (yij)*(H[i]*Hss[j] + Hs[i]*Hs[j]);

				A[0][1] -= (yij)*(Hrs[i]*H[j] + Hr[i]*Hs[j]);
				A[1][0] -= (yij)*(Hrs[i]*H[j] + Hs[i]*Hr[j]);
			}
		}
	
		// determinant of A
		D = A[0][0]*A[1][1] - A[0][1]*A[1][0];

		// solve for u = A^(-1)*R
		u[0] = (A[1][1]*Q[0] - A[0][1]*Q[1])/D;
		u[1] = (A[0][0]*Q[1] - A[1][0]*Q[0])/D;

		// calculate displacement norm
		normu = u[0]*u[0]+u[1]*u[1];

		// check for convergence
		bconv = ((normu < 1e-10));
		if (!bconv && (n <= NMAX))
		{
			// Don't update if converged otherwise the point q
			// does not correspond with the current values for (r,s)
			r += u[0];
			s += u[1];
			++n;
		}
		else break;
	}
	while (1);

	// evaluate q
	q = vec3d(0,0,0);
	for (int i=0; i<ne; ++i) q += y[i]*H[i];

	return bconv;
}

//-----------------------------------------------------------------------------
vec3d FESurface::Position(FESurfaceElement& el, double r, double s)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pMesh;

	// number of element nodes
	int ne = el.Nodes();

	// get the elements nodal positions
	vec3d y[FEElement::MAX_NODES];
    if (!m_bshellb) for (int i = 0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;
    else for (int i = 0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).st();

	double H[FEElement::MAX_NODES];
	el.shape_fnc(H, r, s);

	vec3d q(0,0,0);
	for (int i=0; i<ne; ++i)
	{
		q += y[i]*H[i];
	}

	return q;
}

//-----------------------------------------------------------------------------
//! This function calculates the position of integration point n

vec3d FESurface::Position(FESurfaceElement &el, int n)
{
    // get the mesh to which this surface belongs
    FEMesh& mesh = *m_pMesh;
    
    // number of element nodes
    int ne = el.Nodes();
    
    // get the elements nodal positions
    vec3d y[FEElement::MAX_NODES];
    if (!m_bshellb) for (int i = 0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;
    else for (int i = 0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).st();
    
    double* H = el.H(n);
    
    vec3d q(0,0,0);
    for (int i=0; i<ne; ++i)
    {
        q += y[i]*H[i];
    }
    
    return q;
}

//-----------------------------------------------------------------------------
void FESurface::NodalCoordinates(FESurfaceElement& el, vec3d* re)
{
	int ne = el.Nodes();
	if (!m_bshellb) for (int i = 0; i < ne; ++i) re[i] = Node(el.m_lnode[i]).m_rt;
    else for (int i = 0; i < ne; ++i) re[i] = Node(el.m_lnode[i]).st();
}

//-----------------------------------------------------------------------------
//! detect face normal relative to element:
//! return +1 if face points away from element, -1 if face points into element, 0 if invalid solution found
double FESurface::FacePointing(FESurfaceElement& se, FEElement& el)
{
    FEMesh& mesh = *GetMesh();
    // get point on surface element
    vec3d sp = Position(se, 0,0);
    
    // get surface normal at that point;
    vec3d sn = SurfaceNormal(se, 0,0);
    
    // check if element attached to this surface element is solid or shell
    FESolidElement* sel = dynamic_cast<FESolidElement*>(&el);
    FEShellElement* shl = dynamic_cast<FEShellElement*>(&el);
    
    // get centroid of element el
    vec3d c(0,0,0);
    if (sel) {
        for (int i=0; i<sel->Nodes(); ++i) {
            FENode& node = mesh.Node(sel->m_node[i]);
            c += node.m_rt;
        }
        c /= sel->Nodes();
    }
    else if (shl) {
        for (int i=0; i<sel->Nodes(); ++i) {
            FENode& node = mesh.Node(sel->m_node[i]);
            c += node.m_rt;
            c += node.st();
        }
        c /= (2*sel->Nodes());
    }
    else
        return 0;
    
    // project vector from centroid to surface point onto surface normal
    double d = (sp - c)*sn;
    if (d > 0) return 1.0;
    else if (d < 0) return -1.0;
    else return 0.0;
}

//-----------------------------------------------------------------------------
//! This function calculates the projection of x on the surface element el.
//! It does this by finding the solution of the nonlinear equation (x-y)*y,[a]=0,
//! where the comma denotes differentation and a ranges from 1 to 2.
//! The system is solved using the Newton-Raphson method.
//! The surface element may be either a quad or a triangular element.

vec3d FESurface::ProjectToSurface(FESurfaceElement& el, vec3d x, double& r, double& s)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pMesh;

	// number of element nodes
	int ne = el.Nodes();

	// get the elements nodal positions
	vec3d y[FEElement::MAX_NODES];
    if (!m_bshellb) for (int i=0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;
    else for (int i=0; i<ne; ++i) y[i] = mesh.Node(el.m_node[i]).st();

	// calculate normal projection of x onto element
	vec3d q;
	switch (ne)
	{
	case 3: q = project2tri(y, x, r, s); break;
	case 4: 
		// see if we get lucky
		if (project2quad(y, x, r, s, q) == false)
		{
			// the direct projection failed, so we'll try it more incrementally
			vec3d x0 = (y[0]+y[1]+y[2]+y[3])*0.25;
			r = s = 0;
			bool b = project2quad(y, x0, r, s, q);
			assert(b);
			
			double w = 0.5;
			int l = 1, N = 0, NMAX = 20;
			do
			{
				vec3d xi = x0*(1.0 - w) + x*w;
				b = project2quad(y, xi, r, s, q);
				if (b)
				{
					--l;
					if (l == 0) { x0 = xi; w = 1.0; }
					else w *= 2.0;
				}
				else 
				{
					++l;
					w *= 0.5;
				}
				++N;
			}
			while ((l >= 0) && (N<=NMAX) && (w>0.1));
		}
		break;
	case 6: 
	case 7: 
	case 8:
	case 9:
		if (project2surf(el, y, x, r, s, q)==false)
		{
//			assert(false);
		}
		break;
	default:
		assert(false);
	}

	return q;
}

//-----------------------------------------------------------------------------
//! This function calculates the area of a surface element
double FESurface::FaceArea(FESurfaceElement& el)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pMesh;

	// get the number of nodes
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// get the initial nodes
	vec3d r0[FEElement::MAX_NODES];
	if (!m_bshellb) for (int i=0; i<neln; ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;
    else for (int i=0; i<neln; ++i) r0[i] = mesh.Node(el.m_node[i]).s0();

	// get the integration weights
	double* w = el.GaussWeights();

	double *Gr, *Gs;
	vec3d dxr, dxs;

	double detJ;

	double area = 0;

	int n, k;

	for (n=0; n<nint; ++n)
	{
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		// calculate jacobian
		dxr = dxs = vec3d(0,0,0);
		for (k=0; k<neln; ++k)
		{
			dxr.x += Gr[k]*r0[k].x;
			dxr.y += Gr[k]*r0[k].y;
			dxr.z += Gr[k]*r0[k].z;

			dxs.x += Gs[k]*r0[k].x;
			dxs.y += Gs[k]*r0[k].y;
			dxs.z += Gs[k]*r0[k].z;
		}

		detJ = (dxr ^ dxs).norm();

		area += w[n]*detJ;
	}

	return area;
}

//-----------------------------------------------------------------------------
//! This function calculates the area of a surface element
double FESurface::CurrentFaceArea(FESurfaceElement& el)
{
	// get the mesh to which this surface belongs
	FEMesh& mesh = *m_pMesh;

	// get the number of nodes
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// get the initial nodes
	vec3d rt[FEElement::MAX_NODES];
	if (!m_bshellb) for (int i = 0; i < neln; ++i) rt[i] = mesh.Node(el.m_node[i]).m_rt;
	else for (int i = 0; i < neln; ++i) rt[i] = mesh.Node(el.m_node[i]).st();

	// get the integration weights
	double* w = el.GaussWeights();

	double* Gr, * Gs;
	vec3d dxr, dxs;

	double detJ;

	double area = 0;

	int n, k;

	for (n = 0; n < nint; ++n)
	{
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		// calculate jacobian
		dxr = dxs = vec3d(0, 0, 0);
		for (k = 0; k < neln; ++k)
		{
			dxr.x += Gr[k] * rt[k].x;
			dxr.y += Gr[k] * rt[k].y;
			dxr.z += Gr[k] * rt[k].z;

			dxs.x += Gs[k] * rt[k].x;
			dxs.y += Gs[k] * rt[k].y;
			dxs.z += Gs[k] * rt[k].z;
		}

		detJ = (dxr ^ dxs).norm();

		area += w[n] * detJ;
	}

	return area;
}

//-----------------------------------------------------------------------------
//! Calculate the max element size.
double FESurface::MaxElementSize()
{
	FEMesh& mesh = *m_pMesh;
	double h2 = 0.0;
	int NS = Elements();
	for (int m=0; m<NS; ++m)
	{
		FESurfaceElement& e = Element(m);
		int n = e.Nodes();
		for (int i=0; i<n; ++i)
			for (int j=i+1; j<n; ++j)
			{
                vec3d& a = mesh.Node(e.m_node[i]).m_rt;
                vec3d& b = mesh.Node(e.m_node[j]).m_rt;
				double L2 = (b - a)*(b - a);
				if (L2 > h2) h2 = L2;
			}
	}
	return sqrt(h2);
}

//-----------------------------------------------------------------------------
//! Calculates the metric tensor at the point with surface coordinates (r,s)
//! \todo Perhaps I should place this function in the element class.

mat2d FESurface::Metric0(FESurfaceElement& el, double r, double s)
{
	// nr of element nodes
	int neln = el.Nodes();
	
	// element nodes
	vec3d r0[FEElement::MAX_NODES];
	if (!m_bshellb) for (int i=0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
    else for (int i=0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).s0();
	
	// shape function derivatives
	double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
	
	// get the shape function values at this node
	el.shape_deriv(Hr, Hs, r, s);
	
	// get the tangent vectors
	vec3d t1(0,0,0);
	vec3d t2(0,0,0);
	for (int k=0; k<neln; ++k)
	{
		t1 += r0[k]*Hr[k];
		t2 += r0[k]*Hs[k];
	}
	
	// calculate metric tensor
	return mat2d(t1*t1, t1*t2, t2*t1, t2*t2);
}

//-----------------------------------------------------------------------------
//! Calculates the metric tensor at the point with surface coordinates (r,s)
//! \todo Perhaps I should place this function in the element class.

mat2d FESurface::Metric(FESurfaceElement& el, double r, double s)
{
	// nr of element nodes
	int neln = el.Nodes();
	
	// element nodes
	vec3d rt[FEElement::MAX_NODES];
    if (!m_bshellb) for (int i=0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
    else for (int i=0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).st();
	
	// shape function derivatives
	double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
	
	// get the shape function values at this node
	el.shape_deriv(Hr, Hs, r, s);
	
	// get the tangent vectors
	vec3d t1(0,0,0);
	vec3d t2(0,0,0);
	for (int k=0; k<neln; ++k)
	{
		t1 += rt[k]*Hr[k];
		t2 += rt[k]*Hs[k];
	}
	
	// calculate metric tensor
	return mat2d(t1*t1, t1*t2, t2*t1, t2*t2);
}

//-----------------------------------------------------------------------------
//! Calculates the metric tensor at an integration point
//! \todo Perhaps I should place this function in the element class.

mat2d FESurface::Metric(const FESurfaceElement& el, int n) const
{
    // nr of element nodes
    int neln = el.Nodes();
    
    // element nodes
    vec3d rt[FEElement::MAX_NODES];
    if (!m_bshellb) for (int i=0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
    else for (int i=0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).st();
    
    // get the shape function derivatives at this integration point
    double* Hr = el.Gr(n);
    double* Hs = el.Gs(n);
    
    // get the tangent vectors
    vec3d t1(0,0,0);
    vec3d t2(0,0,0);
    for (int k=0; k<neln; ++k)
    {
        t1 += rt[k]*Hr[k];
        t2 += rt[k]*Hs[k];
    }
    
    // calculate metric tensor
    return mat2d(t1*t1, t1*t2, t2*t1, t2*t2);
}

//-----------------------------------------------------------------------------
//! Calculates the metric tensor at an integration point at previous time
//! \todo Perhaps I should place this function in the element class.

mat2d FESurface::MetricP(FESurfaceElement& el, int n)
{
    // nr of element nodes
    int neln = el.Nodes();
    
    // element nodes
    vec3d rp[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i) rp[i] = m_pMesh->Node(el.m_node[i]).m_rp;
    
    // get the shape function derivatives at this integration point
    double* Hr = el.Gr(n);
    double* Hs = el.Gs(n);
    
    // get the tangent vectors
    vec3d t1(0,0,0);
    vec3d t2(0,0,0);
    for (int k=0; k<neln; ++k)
    {
        t1 += rp[k]*Hr[k];
        t2 += rp[k]*Hs[k];
    }
    
    // calculate metric tensor
    return mat2d(t1*t1, t1*t2, t2*t1, t2*t2);
}

//-----------------------------------------------------------------------------
//! Given an element an the natural coordinates of a point in this element, this
//! function returns the global position vector.
vec3d FESurface::Local2Global(FESurfaceElement &el, double r, double s)
{
	// get the mesh
	FEMesh& mesh = *m_pMesh;

	// get the coordinates of the element nodes
	int ne = el.Nodes();
	vec3d y[FEElement::MAX_NODES];
    if (!m_bshellb) for (int l=0; l<ne; ++l) y[l] = mesh.Node(el.m_node[l]).m_rt;
    else for (int l=0; l<ne; ++l) y[l] = mesh.Node(el.m_node[l]).st();

	// calculate the element position
	return el.eval(y, r, s);
}

//-----------------------------------------------------------------------------
//! This function calculates the global location of an integration point
//!

vec3d FESurface::Local2Global(FESurfaceElement &el, int n)
{
	FEMesh& m = *m_pMesh;

	// get the shape functions at this integration point
	double* H = el.H(n);

	// calculate the location
	vec3d r(0);
	int ne = el.Nodes();
    if (!m_bshellb) for (int i=0; i<ne; ++i) r += m.Node(el.m_node[i]).m_rt*H[i];
    else for (int i=0; i<ne; ++i) r += m.Node(el.m_node[i]).st()*H[i];

	return r;
}

//-----------------------------------------------------------------------------
//! Given an element an the natural coordinates of a point in this element, this
//! function returns the global position vector at the previous time.
vec3d FESurface::Local2GlobalP(FESurfaceElement &el, double r, double s)
{
    // get the mesh
    FEMesh& mesh = *m_pMesh;
    
    // get the coordinates of the element nodes
    int ne = el.Nodes();
    vec3d y[FEElement::MAX_NODES];
    for (int l=0; l<ne; ++l) y[l] = mesh.Node(el.m_node[l]).m_rp;
    
    // calculate the element position
    return el.eval(y, r, s);
}

//-----------------------------------------------------------------------------
//! This function calculates the global location of an integration point
//!

vec3d FESurface::Local2GlobalP(FESurfaceElement &el, int n)
{
    FEMesh& m = *m_pMesh;
    
    // get the shape functions at this integration point
    double* H = el.H(n);
    
    // calculate the location
    vec3d r;
    int ne = el.Nodes();
    for (int i=0; i<ne; ++i) r += m.Node(el.m_node[i]).m_rp*H[i];
    
    return r;
}

//-----------------------------------------------------------------------------
//! This function calculates the normal of a surface element at integration
//! point n

vec3d FESurface::SurfaceNormal(const FESurfaceElement &el, int n) const
{
	FEMesh& m = *m_pMesh;

	// get the shape function derivatives at this integration point
	double* Hr = el.Gr(n);
	double* Hs = el.Gs(n);

	// get the coordinates of the element nodes
	int ne = el.Nodes();
	vec3d y[FEElement::MAX_NODES];
    if (!m_bshellb) for (int i=0; i<ne; ++i) y[i] = m.Node(el.m_node[i]).m_rt;
    else for (int i=0; i<ne; ++i) y[i] = m.Node(el.m_node[i]).st();

	// calculate the tangents
	vec3d xr, xs;
	for (int i=0; i<ne; ++i)
	{
		xr += y[i]*Hr[i];
		xs += y[i]*Hs[i];
	}

	// calculate the normal
	vec3d np = xr ^ xs;
	np.unit();
    if (m_bshellb) np = -np;

	return np;
}

//-----------------------------------------------------------------------------
//! This function calculates the normal of a surface element at the natural
//! coordinates (r,s)

vec3d FESurface::SurfaceNormal(FESurfaceElement &el, double r, double s) const
{
	int l;
	FEMesh& mesh = *m_pMesh;
	
	// get the coordinates of the element nodes
	int ne = el.Nodes();
	vec3d y[FEElement::MAX_NODES];
    if (!m_bshellb) for (l=0; l<ne; ++l) y[l] = mesh.Node(el.m_node[l]).m_rt;
    else for (l=0; l<ne; ++l) y[l] = mesh.Node(el.m_node[l]).st();
	
	// set up shape functions and derivatives
	double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
	el.shape_deriv(Hr, Hs, r, s);
	
	// calculate the element tangents
	vec3d xr(0,0,0), xs;
	for (l=0; l<ne; ++l)
	{
		xr += y[l]*Hr[l];
		xs += y[l]*Hs[l];
	}
	
	// calculate the normal
	vec3d np = xr ^ xs;
	np.unit();
    
    if (m_bshellb) np = -np;
	
	return np;
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the
//! element normals at the node

void FESurface::UpdateNodeNormals()
{
    const int MN = FEElement::MAX_NODES;
    vec3d y[MN];
    
    // zero nodal normals
    zero(m_nn);
    
    // loop over all elements
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int ne = el.Nodes();
        
        // get the nodal coordinates
        for (int j=0; j<ne; ++j) y[j] = Node(el.m_lnode[j]).m_rt;
        
        // calculate the normals
        for (int j=0; j<ne; ++j)
        {
            int jp1 = (j+1)%ne;
            int jm1 = (j+ne-1)%ne;
            vec3d n = (y[jp1] - y[j]) ^ (y[jm1] - y[j]);
            m_nn[el.m_lnode[j]] += n;
        }
    }
    
    // normalize all vectors
    const int N = Nodes();
    for (int i=0; i<N; ++i) m_nn[i].unit();
}

//-----------------------------------------------------------------------------
//! This function checks whether the point with natural coordinates (r,s) is
//! inside the element, within a tolerance of tol.

bool FESurface::IsInsideElement(FESurfaceElement& el, double r, double s, double tol)
{
	int ne = el.Nodes();
	switch (ne)
	{
	case 4:
	case 8:
	case 9:
		{
			// check quads
			if ((r >= -1-tol) && (r <= 1+tol) && (s >= -1-tol) && (s <= 1+tol)) return true;
			else return false;
		}
		break;
	case 3:
	case 6:
	case 7:
		{
			// check triangles
			if ((r >= -tol) && (s >= -tol) && (r+s <= 1+tol)) return true;
			else return false;
		}
		break;
	}
	assert(false);
	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant base vectors of a surface element
//! at an integration point

void FESurface::CoBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2])
{
	FEMesh& m = *m_pMesh;

	// get the nr of nodes
	int n = el.Nodes();

	// get the shape function derivatives
	double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
	el.shape_deriv(Hr, Hs, r, s);

	t[0] = t[1] = vec3d(0,0,0);
    if (!m_bshellb) {
        for (int i=0; i<n; ++i)
        {
            t[0] += m.Node(el.m_node[i]).m_rt*Hr[i];
            t[1] += m.Node(el.m_node[i]).m_rt*Hs[i];
        }
    }
    else {
        for (int i=0; i<n; ++i)
        {
            t[0] -= m.Node(el.m_node[i]).st()*Hr[i];
            t[1] -= m.Node(el.m_node[i]).st()*Hs[i];
        }
    }
}


//-----------------------------------------------------------------------------
//! This function calculates the covariant base vectors of a surface element
//! at an integration point

void FESurface::CoBaseVectors(const FESurfaceElement& el, int j, vec3d t[2]) const
{
	FEMesh& m = *m_pMesh;

	// get the nr of nodes
	int n = el.Nodes();

	// get the shape function derivatives
	double* Hr = el.Gr(j);
	double* Hs = el.Gs(j);

	t[0] = t[1] = vec3d(0,0,0);
    if (!m_bshellb) {
        for (int i=0; i<n; ++i)
        {
            t[0] += m.Node(el.m_node[i]).m_rt*Hr[i];
            t[1] += m.Node(el.m_node[i]).m_rt*Hs[i];
        }
    }
    else {
        for (int i=0; i<n; ++i)
        {
            t[0] -= m.Node(el.m_node[i]).st()*Hr[i];
            t[1] -= m.Node(el.m_node[i]).st()*Hs[i];
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant base vectors of a surface element
//! at an integration point at previous time step

void FESurface::CoBaseVectorsP(FESurfaceElement& el, int j, vec3d t[2])
{
    FEMesh& m = *m_pMesh;
    
    // get the nr of nodes
    int n = el.Nodes();
    
    // get the shape function derivatives
    double* Hr = el.Gr(j);
    double* Hs = el.Gs(j);
    
    t[0] = t[1] = vec3d(0,0,0);
    for (int i=0; i<n; ++i)
    {
        t[0] += m.Node(el.m_node[i]).m_rp*Hr[i];
        t[1] += m.Node(el.m_node[i]).m_rp*Hs[i];
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant base vectors of a surface element
//! at the natural coordinates (r,s)

void FESurface::CoBaseVectors0(FESurfaceElement &el, double r, double s, vec3d t[2])
{
	int i;
	const int MN = FEElement::MAX_NODES;
	vec3d y[MN];
	double H0[MN], H1[MN];
	int n = el.Nodes();
	if (!m_bshellb) for (i=0; i<n; ++i) y[i] = m_pMesh->Node(el.m_node[i]).m_r0;
    else for (i=0; i<n; ++i) y[i] = m_pMesh->Node(el.m_node[i]).s0();
	el.shape_deriv(H0, H1, r, s);
	t[0] = t[1] = vec3d(0,0,0);
	for (i=0; i<n; ++i) 
	{
		t[0] += y[i]*H0[i];
		t[1] += y[i]*H1[i];
	}
}

//-----------------------------------------------------------------------------
void FESurface::ContraBaseVectors(const FESurfaceElement& el, int j, vec3d t[2]) const
{
    vec3d e[2];
    CoBaseVectors(el, j, e);
    mat2d M = Metric(el, j);
    mat2d Mi = M.inverse();
    
    t[0] = e[0]*Mi[0][0] + e[1]*Mi[0][1];
    t[1] = e[0]*Mi[1][0] + e[1]*Mi[1][1];
}

//-----------------------------------------------------------------------------
void FESurface::ContraBaseVectorsP(FESurfaceElement& el, int j, vec3d t[2])
{
    vec3d e[2];
    CoBaseVectorsP(el, j, e);
    mat2d M = MetricP(el, j);
    mat2d Mi = M.inverse();
    
    t[0] = e[0]*Mi[0][0] + e[1]*Mi[0][1];
    t[1] = e[0]*Mi[1][0] + e[1]*Mi[1][1];
}

//-----------------------------------------------------------------------------
void FESurface::ContraBaseVectors(FESurfaceElement& el, double r, double s, vec3d t[2])
{
	vec3d e[2];
	CoBaseVectors(el, r, s, e);
	mat2d M = Metric(el, r, s);
	mat2d Mi = M.inverse();

	t[0] = e[0]*Mi[0][0] + e[1]*Mi[0][1];
	t[1] = e[0]*Mi[1][0] + e[1]*Mi[1][1];
}

//-----------------------------------------------------------------------------

void FESurface::ContraBaseVectors0(FESurfaceElement& el, double r, double s, vec3d t[2])
{
	vec3d e[2];
	CoBaseVectors0(el, r, s, e);
	mat2d M = Metric0(el, r, s);
	mat2d Mi = M.inverse();

	t[0] = e[0]*Mi[0][0] + e[1]*Mi[0][1];
	t[1] = e[0]*Mi[1][0] + e[1]*Mi[1][1];
}

//-----------------------------------------------------------------------------
double FESurface::jac0(FESurfaceElement &el, int n)
{
	const int nseln = el.Nodes();
	vec3d r0[FEElement::MAX_NODES];
	if (!m_bshellb) for (int i=0; i<nseln; ++i) r0[i] = GetMesh()->Node(el.m_node[i]).m_r0;
    else for (int i=0; i<nseln; ++i) r0[i] = GetMesh()->Node(el.m_node[i]).s0();

	double* Gr = el.Gr(n);
	double* Gs = el.Gs(n);
	vec3d dxr(0,0,0), dxs(0,0,0);
	for (int k=0; k<nseln; ++k)
	{
		dxr.x += Gr[k]*r0[k].x;
		dxr.y += Gr[k]*r0[k].y;
		dxr.z += Gr[k]*r0[k].z;

		dxs.x += Gs[k]*r0[k].x;
		dxs.y += Gs[k]*r0[k].y;
		dxs.z += Gs[k]*r0[k].z;
	}
	return (dxr ^ dxs).norm();
}

//-----------------------------------------------------------------------------
double FESurface::jac0(const FESurfaceElement &el, int n, vec3d& nu)
{
	const int nseln = el.Nodes();
	vec3d r0[FEElement::MAX_NODES];
	if (!m_bshellb) for (int i=0; i<nseln; ++i) r0[i] = GetMesh()->Node(el.m_node[i]).m_r0;
    else for (int i=0; i<nseln; ++i) r0[i] = GetMesh()->Node(el.m_node[i]).s0();

	double* Gr = el.Gr(n);
	double* Gs = el.Gs(n);
	vec3d dxr(0,0,0), dxs(0,0,0);
	for (int k=0; k<nseln; ++k)
	{
		dxr.x += Gr[k]*r0[k].x;
		dxr.y += Gr[k]*r0[k].y;
		dxr.z += Gr[k]*r0[k].z;

		dxs.x += Gs[k]*r0[k].x;
		dxs.y += Gs[k]*r0[k].y;
		dxs.z += Gs[k]*r0[k].z;
	}
	nu = dxr ^ dxs;
	return nu.unit();
}

//-----------------------------------------------------------------------------
// This function calculates the intersection of a ray with a triangle
// and returns true if the ray intersects the triangle.
//
bool IntersectTri(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	vec3d e[2], E[2];
	mat2d G;

	// create base vectors on triangle
	e[0] = y[1]-y[0];
	e[1] = y[2]-y[0];

	// create triangle normal
	vec3d m = e[0]^e[1]; m.unit();

	double d = n*m;
//	if (d != 0)
	// normals should be pointing opposite to each other for valid contact
	if (d < 0)
	{
		// distance from r to plane of triangle
		g = m*(y[0] - r)/d;

		// intersection point with plane of triangle
		vec3d q = r + n*g;

		// next, we decompose q into its components
		// in the triangle basis
		// we need to create the dual basis
		// first, we calculate the metric tensor
		G[0][0] = e[0]*e[0]; G[0][1] = e[0]*e[1];
		G[1][0] = e[1]*e[0]; G[1][1] = e[1]*e[1];

		// and its inverse
		mat2d Gi = G.inverse();

		// build dual basis
		E[0] = e[0]*Gi[0][0] + e[1]*Gi[0][1];
		E[1] = e[0]*Gi[1][0] + e[1]*Gi[1][1];

		// get the components
		rs[0] = E[0]*(q - y[0]);
		rs[1] = E[1]*(q - y[0]);

		// see if the intersection point is inside the triangle
		if ((rs[0] >= -eps) && (rs[1] >= -eps) && (rs[0]+rs[1] <= 1+eps)) return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a quad
//! and returns true if the ray intersected.
//!
bool IntersectQuad(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	// first we're going to see if the ray intersects the two subtriangles
	vec3d x1[3], x2[3];
	x1[0] = y[0]; x2[0] = y[2];
	x1[1] = y[1]; x2[1] = y[3];
	x1[2] = y[3]; x2[2] = y[1];

	bool b = false;
	double rp = 0, sp = 0;

	if (IntersectTri(x1, r, n, rs, g, eps))
	{
		// we've intersected the first triangle
		b = true;
		rp = -1.0 + 2.0*rs[0];
		sp = -1.0 + 2.0*rs[1];
	}
	else if (IntersectTri(x2, r, n, rs, g, eps))
	{
		// we've intersected the second triangle
		b = true;
		rp = 1.0 - 2.0*rs[0];
		sp = 1.0 - 2.0*rs[1];
	}

	// if one of the triangels was intersected,
	// we calculate a more accurate projection
	if (b)
	{
		mat3d A;
		vec3d dx;
		vec3d F, F1, F2, F3;
		double H[4], H1[4], H2[4];

		double l1 = rp;
		double l2 = sp;
		double l3 = g;
		
		int nn = 0;
		int maxn = 5;
		do
		{
			// shape functions of quad
			H[0] = 0.25*(1 - l1)*(1 - l2);
			H[1] = 0.25*(1 + l1)*(1 - l2);
			H[2] = 0.25*(1 + l1)*(1 + l2);
			H[3] = 0.25*(1 - l1)*(1 + l2);

			// shape function derivatives
			H1[0] = -0.25*(1 - l2); H2[0] = -0.25*(1 - l1);
			H1[1] =  0.25*(1 - l2); H2[1] = -0.25*(1 + l1);
			H1[2] =  0.25*(1 + l2); H2[2] =  0.25*(1 + l1);
			H1[3] = -0.25*(1 + l2); H2[3] =  0.25*(1 - l1);

			// calculate residual
			F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3];

			// residual derivatives
			F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3];
			F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3];
			F3 = n;

			// set up the tangent matrix
			A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
			A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
			A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;

			// calculate solution increment
            mat3d Ai = A.inverse();
			dx = -(Ai*F);

			// update solution
			l1 += dx.x;
			l2 += dx.y;
			l3 += dx.z;

			++nn;
		}
		while ((dx.norm() > 1e-7) && (nn < maxn));
        if (nn == maxn) return false;

		// store results
		rs[0] = l1;
		rs[1] = l2;
		g     = l3;
        vec3d nu2 = F1 ^ F2;
        nu2.unit();
        double cosq = n*nu2;
        if (cosq > 0) return false;

		// see if the point is inside the quad
		if ((rs[0] >= -1-eps) && (rs[0] <= 1+eps) && 
			(rs[1] >= -1-eps) && (rs[1] <= 1+eps)) return true;
	}

	return false;
}


//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a quad
//! and returns true if the ray intersected.
//!
bool IntersectQuad8(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	mat3d A;
	vec3d dx;
	vec3d F, F1, F2, F3;
	double H[8], H1[8], H2[8];

	double l1 = 0;
	double l2 = 0;
	double l3 = 0;
		
	int nn = 0;
	int maxn = 10;
	do
	{
		// shape functions of quad
		H[4] = 0.5*(1 - l1*l1)*(1 - l2);
		H[5] = 0.5*(1 - l2*l2)*(1 + l1);
		H[6] = 0.5*(1 - l1*l1)*(1 + l2);
		H[7] = 0.5*(1 - l2*l2)*(1 - l1);

		H[0] = 0.25*(1 - l1)*(1 - l2) - 0.5*(H[4] + H[7]);
		H[1] = 0.25*(1 + l1)*(1 - l2) - 0.5*(H[4] + H[5]);
		H[2] = 0.25*(1 + l1)*(1 + l2) - 0.5*(H[5] + H[6]);
		H[3] = 0.25*(1 - l1)*(1 + l2) - 0.5*(H[6] + H[7]);

		// shape function derivatives
		H1[4] = -l1*(1 - l2);
		H1[5] = 0.5*(1 - l2*l2);
		H1[6] = -l1*(1 + l2);
		H1[7] = -0.5*(1 - l2*l2);

		H1[0] = -0.25*(1 - l2) - 0.5*(H1[4] + H1[7]);
		H1[1] =  0.25*(1 - l2) - 0.5*(H1[4] + H1[5]);
		H1[2] =  0.25*(1 + l2) - 0.5*(H1[5] + H1[6]);
		H1[3] = -0.25*(1 + l2) - 0.5*(H1[6] + H1[7]);

		H2[4] = -0.5*(1 - l1*l1);
		H2[5] = -l2*(1 + l1);
		H2[6] = 0.5*(1 - l1*l1);
		H2[7] = -l2*(1 - l1);

		H2[0] = -0.25*(1 - l1) - 0.5*(H2[4] + H2[7]);
		H2[1] = -0.25*(1 + l1) - 0.5*(H2[4] + H2[5]);
		H2[2] =  0.25*(1 + l1) - 0.5*(H2[5] + H2[6]);
		H2[3] =  0.25*(1 - l1) - 0.5*(H2[6] + H2[7]);
		
		// calculate residual
		F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3] - y[4]*H[4] - y[5]*H[5] - y[6]*H[6] - y[7]*H[7];

		// residual derivatives
		F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3] - y[4]*H1[4] - y[5]*H1[5] - y[6]*H1[6] - y[7]*H1[7];
		F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3] - y[4]*H2[4] - y[5]*H2[5] - y[6]*H2[6] - y[7]*H2[7];
		F3 = n;

		// set up the tangent matrix
		A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
		A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
		A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;

		// calculate solution increment
        mat3d Ai = A.inverse();
		dx = -(Ai*F);

		// update solution
		l1 += dx.x;
		l2 += dx.y;
		l3 += dx.z;

		++nn;
	}
	while ((dx.norm() > 1e-7) && (nn < maxn));
    if (nn == maxn) return false;

	// store results
	rs[0] = l1;
	rs[1] = l2;
	g     = l3;
    vec3d nu2 = F1 ^ F2;
    nu2.unit();
    double cosq = n*nu2;
    if (cosq > 0) return false;

	// see if the point is inside the quad
	if ((rs[0] >= -1-eps) && (rs[0] <= 1+eps) && 
		(rs[1] >= -1-eps) && (rs[1] <= 1+eps)) return true;

	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a quad
//! and returns true if the ray intersected.
//!
bool IntersectQuad9(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	mat3d A;
	vec3d dx;
	vec3d F, F1, F2, F3;
	double H[9], H1[9], H2[9];

	double l1 = 0;
	double l2 = 0;
	double l3 = 0;
		
	int nn = 0;
	int maxn = 10;
	do
	{
		// shape functions of quad
		double R[3] = {0.5*l1*(l1-1.0), 0.5*l1*(l1+1.0), 1.0 - l1*l1};
		double S[3] = {0.5*l2*(l2-1.0), 0.5*l2*(l2+1.0), 1.0 - l2*l2};
		double DR[3] = {l1-0.5, l1+0.5, -2.0*l1};
		double DS[3] = {l2-0.5, l2+0.5, -2.0*l2};

		H[0] = R[0]*S[0];
		H[1] = R[1]*S[0];
		H[2] = R[1]*S[1];
		H[3] = R[0]*S[1];
		H[4] = R[2]*S[0];
		H[5] = R[1]*S[2];
		H[6] = R[2]*S[1];
		H[7] = R[0]*S[2];
		H[8] = R[2]*S[2];

		// shape function derivatives
		H1[0] = DR[0]*S[0];
		H1[1] = DR[1]*S[0];
		H1[2] = DR[1]*S[1];
		H1[3] = DR[0]*S[1];
		H1[4] = DR[2]*S[0];
		H1[5] = DR[1]*S[2];
		H1[6] = DR[2]*S[1];
		H1[7] = DR[0]*S[2];
		H1[8] = DR[2]*S[2];

		H2[0] = R[0]*DS[0];
		H2[1] = R[1]*DS[0];
		H2[2] = R[1]*DS[1];
		H2[3] = R[0]*DS[1];
		H2[4] = R[2]*DS[0];
		H2[5] = R[1]*DS[2];
		H2[6] = R[2]*DS[1];
		H2[7] = R[0]*DS[2];
		H2[8] = R[2]*DS[2];
		
		// calculate residual
		F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3] - y[4]*H[4] - y[5]*H[5] - y[6]*H[6] - y[7]*H[7] - y[8]*H[8];

		// residual derivatives
		F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3] - y[4]*H1[4] - y[5]*H1[5] - y[6]*H1[6] - y[7]*H1[7] - y[8]*H1[8];
		F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3] - y[4]*H2[4] - y[5]*H2[5] - y[6]*H2[6] - y[7]*H2[7] - y[8]*H2[8];
		F3 = n;

		// set up the tangent matrix
		A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
		A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
		A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;

		// calculate solution increment
        mat3d Ai = A.inverse();
		dx = -(Ai*F);

		// update solution
		l1 += dx.x;
		l2 += dx.y;
		l3 += dx.z;

		++nn;
	}
	while ((dx.norm() > 1e-7) && (nn < maxn));
    if (nn == maxn) return false;

	// store results
	rs[0] = l1;
	rs[1] = l2;
	g     = l3;
    vec3d nu2 = F1 ^ F2;
    nu2.unit();
    double cosq = n*nu2;
    if (cosq > 0) return false;

	// see if the point is inside the quad
	if ((rs[0] >= -1-eps) && (rs[0] <= 1+eps) && 
		(rs[1] >= -1-eps) && (rs[1] <= 1+eps)) return true;

	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a 6-node triangle
//! and returns true if the ray intersected.
//!
bool IntersectTri6(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	// first we're going to see if the ray intersects the four subtriangles
	vec3d x1[3], x2[3], x3[3], x4[3];
	x1[0] = y[0]; x2[0] = y[3]; x3[0] = y[5]; x4[0] = y[4];
	x1[1] = y[3]; x2[1] = y[1]; x3[1] = y[4]; x4[1] = y[5];
	x1[2] = y[5]; x2[2] = y[4]; x3[2] = y[2]; x4[2] = y[3];
	
	bool b = false;
	double rp = 0, sp = 0;
	
	if (IntersectTri(x1, r, n, rs, g, eps))
	{
		// we've intersected the first triangle
		b = true;
		rp = rs[0]/2.0;
		sp = rs[1]/2.0;
	}
	else if (IntersectTri(x2, r, n, rs, g, eps))
	{
		// we've intersected the second triangle
		b = true;
		rp = 0.5 + rs[0]/2.0;
		sp = rs[1]/2.0;
	}
	else if (IntersectTri(x3, r, n, rs, g, eps))
	{
		// we've intersected the third triangle
		b = true;
		rp = rs[0]/2.0;
		sp = 0.5 + rs[1]/2.0;
	}
	else if (IntersectTri(x4, r, n, rs, g, eps))
	{
		// we've intersected the fourth triangle
		b = true;
		rp = 0.5 - rs[0]/2.0;
		sp = 0.5 - rs[1]/2.0;
	}
	
	// if one of the triangels was intersected,
	// we calculate a more accurate projection
	if (b)
	{
		mat3d A;
		vec3d dx;
		vec3d F, F1, F2, F3;
		double H[6], H1[6], H2[6];
		
		double l0;
		double l1 = rp;
		double l2 = sp;
		double l3 = g;
		
		int nn = 0;
		int maxn = 10;
		do
		{
			l0 = 1 - l1 - l2;
			
			// shape functions of quad
			H[0] = l0*(2.0*l0 - 1.0);
			H[1] = l1*(2.0*l1 - 1.0);
			H[2] = l2*(2.0*l2 - 1.0);
			H[3] = 4.0*l0*l1;
			H[4] = 4.0*l1*l2;
			H[5] = 4.0*l2*l0;
			
			// shape function derivatives
			H1[0] = -3.0 + 4.0*l1 + 4.0*l2;
			H1[1] =  4.0*l1 - 1.0;
			H1[2] =  0.0;
			H1[3] =  4.0 - 8.0*l1 - 4.0*l2;
			H1[4] =  4.0*l2;
			H1[5] = -4.0*l2;
			
			H2[0] = -3.0 + 4.0*l2 + 4.0*l1;
			H2[1] =  0.0;
			H2[2] =  4.0*l2 - 1.0;
			H2[3] = -4.0*l1;
			H2[4] =  4.0*l1;
			H2[5] =  4.0 - 8.0*l2 - 4.0*l1;

			// calculate residual
			F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3] - y[4]*H[4] - y[5]*H[5];
			
			// residual derivatives
			F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3] - y[4]*H1[4] - y[5]*H1[5];
			F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3] - y[4]*H2[4] - y[5]*H2[5];
			F3 = n;
			
			// set up the tangent matrix
			A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
			A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
			A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;
			
			// calculate solution increment
            mat3d Ai = A.inverse();
			dx = -(Ai*F);
			
			// update solution
			l1 += dx.x;
			l2 += dx.y;
			l3 += dx.z;
			
			++nn;
		}
		while ((dx.norm() > 1e-7) && (nn < maxn));
        if (nn == maxn) return false;

		// store results
		rs[0] = l1;
		rs[1] = l2;
		g     = l3;
        vec3d nu2 = F1 ^ F2;
        nu2.unit();
        double cosq = n*nu2;
        if (cosq > 0) return false;

		// see if the point is inside the quad
		if ((rs[0] >= -eps) && (rs[1] >= -eps) && (rs[0]+rs[1] <= 1+eps)) return true;
	}
	
	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a 6-node triangle
//! and returns true if the ray intersected.
//!
bool IntersectTri7(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	// first we're going to see if the ray intersects the four subtriangles
	vec3d x1[3], x2[3], x3[3], x4[3];
	x1[0] = y[0]; x2[0] = y[3]; x3[0] = y[5]; x4[0] = y[4];
	x1[1] = y[3]; x2[1] = y[1]; x3[1] = y[4]; x4[1] = y[5];
	x1[2] = y[5]; x2[2] = y[4]; x3[2] = y[2]; x4[2] = y[3];
	
	bool b = false;
	double rp = 0, sp = 0;
	
	if (IntersectTri(x1, r, n, rs, g, eps))
	{
		// we've intersected the first triangle
		b = true;
		rp = rs[0]/2.0;
		sp = rs[1]/2.0;
	}
	else if (IntersectTri(x2, r, n, rs, g, eps))
	{
		// we've intersected the second triangle
		b = true;
		rp = 0.5 + rs[0]/2.0;
		sp = rs[1]/2.0;
	}
	else if (IntersectTri(x3, r, n, rs, g, eps))
	{
		// we've intersected the third triangle
		b = true;
		rp = rs[0]/2.0;
		sp = 0.5 + rs[1]/2.0;
	}
	else if (IntersectTri(x4, r, n, rs, g, eps))
	{
		// we've intersected the fourth triangle
		b = true;
		rp = 0.5 - rs[0]/2.0;
		sp = 0.5 - rs[1]/2.0;
	}
	
	// if one of the triangels was intersected,
	// we calculate a more accurate projection
	if (b)
	{
		mat3d A;
		vec3d dx;
		vec3d F, F1, F2, F3;
		double H[7], H1[7], H2[7];
		
		double l0;
		double l1 = rp;
		double l2 = sp;
		double l3 = g;
		
		int nn = 0;
		int maxn = 5;
		do
		{
			l0 = 1 - l1 - l2;
			
			H[6] = 27.0*l0*l1*l2;
			H[0] = l0*(2.0*l0 - 1.0) + H[6]/9.0;
			H[1] = l1*(2.0*l1 - 1.0) + H[6]/9.0;
			H[2] = l2*(2.0*l2 - 1.0) + H[6]/9.0;
			H[3] = 4.0*l0*l1 - 4.0*H[6]/9.0;
			H[4] = 4.0*l1*l2 - 4.0*H[6]/9.0;
			H[5] = 4.0*l2*l0 - 4.0*H[6]/9.0;

			
			// shape function derivatives
			H1[6] = 27.0*l2*(1.0 - 2.0*l1 - l2);
			H1[0] = -3.0 + 4.0*l1 + 4.0*l2 +     H1[6]/9.0;
			H1[1] =  4.0*l1 - 1.0          +     H1[6]/9.0;
			H1[2] =  0.0                   +     H1[6]/9.0;
			H1[3] =  4.0 - 8.0*l1 - 4.0*l2 - 4.0*H1[6]/9.0;
			H1[4] =  4.0*l2                - 4.0*H1[6]/9.0;
			H1[5] = -4.0*l2                - 4.0*H1[6]/9.0;

			H2[6] = 27.0*l1*(1.0 - l1 - 2.0*l2);
			H2[0] = -3.0 + 4.0*l2 + 4.0*l1     + H2[6]/9.0;
			H2[1] =  0.0                       + H2[6]/9.0;
			H2[2] =  4.0*l2 - 1.0              + H2[6]/9.0;
			H2[3] = -4.0*l1                - 4.0*H2[6]/9.0;
			H2[4] =  4.0*l1                - 4.0*H2[6]/9.0;
			H2[5] =  4.0 - 8.0*l2 - 4.0*l1 - 4.0*H2[6]/9.0;

			// calculate residual
			F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3] - y[4]*H[4] - y[5]*H[5] - y[6]*H[6];
			
			// residual derivatives
			F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3] - y[4]*H1[4] - y[5]*H1[5] - y[6]*H1[6];
			F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3] - y[4]*H2[4] - y[5]*H2[5] - y[6]*H2[6];
			F3 = n;
			
			// set up the tangent matrix
			A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
			A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
			A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;
			
			// calculate solution increment
            mat3d Ai = A.inverse();
			dx = -(Ai*F);
			
			// update solution
			l1 += dx.x;
			l2 += dx.y;
			l3 += dx.z;
			
			++nn;
		}
		while ((dx.norm() > 1e-7) && (nn < maxn));
        if (nn == maxn) return false;

		// store results
		rs[0] = l1;
		rs[1] = l2;
		g     = l3;
        vec3d nu2 = F1 ^ F2;
        nu2.unit();
        double cosq = n*nu2;
        if (cosq > 0) return false;

		// see if the point is inside the quad
		if ((rs[0] >= -eps) && (rs[1] >= -eps) && (rs[0]+rs[1] <= 1+eps)) return true;
	}
	
	return false;
}

//-----------------------------------------------------------------------------
void FESurface::Invert()
{
	ForEachSurfaceElement([](FESurfaceElement& el) {
		int tmp;
		switch (el.Shape())
		{
		case ET_TRI3 : tmp = el.m_node[1]; el.m_node[1] = el.m_node[2]; el.m_node[2] = tmp; break;
		case ET_QUAD4: tmp = el.m_node[1]; el.m_node[1] = el.m_node[3]; el.m_node[3] = tmp; break;
		case ET_QUAD8:
		case ET_QUAD9:
			tmp = el.m_node[1]; el.m_node[1] = el.m_node[3]; el.m_node[3] = tmp;
			tmp = el.m_node[4]; el.m_node[4] = el.m_node[7]; el.m_node[7] = tmp;
			tmp = el.m_node[5]; el.m_node[5] = el.m_node[6]; el.m_node[6] = tmp;
			break;
		case ET_TRI6:
		case ET_TRI7:
			tmp = el.m_node[1]; el.m_node[1] = el.m_node[2]; el.m_node[2] = tmp;
			tmp = el.m_node[3]; el.m_node[3] = el.m_node[5]; el.m_node[5] = tmp;
			break;
		default:
			assert(false);
		}
	});
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a surface element.
//! It simply calls the correct intersection function based on the type
//! of element.
//!
bool FESurface::Intersect(FESurfaceElement& el, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	int N = el.Nodes();

	// get the element nodes
	FEMesh& mesh = *m_pMesh;
	vec3d y[FEElement::MAX_NODES];
    if (!m_bshellb) for (int i=0; i<N; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;
    else for (int i=0; i<N; ++i) y[i] = mesh.Node(el.m_node[i]).st();

	// call the correct intersection function
	switch (N)
	{
	case 3: return IntersectTri  (y, r, n, rs, g, eps); break;
	case 4: return IntersectQuad (y, r, n, rs, g, eps); break;
	case 6: return IntersectTri6 (y, r, n, rs, g, eps); break;
	case 7: return IntersectTri7 (y, r, n, rs, g, eps); break;
	case 8: return IntersectQuad8(y, r, n, rs, g, eps); break;
	case 9: return IntersectQuad9(y, r, n, rs, g, eps); break;
	default:
		assert(false);
	}

	// if we get here, the ray did not intersect the element
	return false;
}

//-----------------------------------------------------------------------------
void FESurface::Serialize(DumpStream &ar)
{
	FEMeshPartition::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		ar & m_surf;
		ar & m_bitfc;
		ar & m_alpha;
		ar & m_bshellb;
		ar & m_el;

		// reallocate integration point data on loading
		if (ar.IsSaving() == false)
		{
			for (int i = 0; i < Elements(); ++i)
			{
				FESurfaceElement& el = Element(i);
				el.SetMeshPartition(this);
				int nint = el.GaussPoints();
				for (int n = 0; n < nint; ++n)
				{
					FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(CreateMaterialPoint());
					assert(pt);
					el.SetMaterialPointData(pt, n);
				}
			}
		}
	}

	if ((ar.IsShallow() == false) && (ar.IsSaving() == false))
	{
		// see if we can find all elements that the faces belong to
		int ne = Elements();
		for (int i = 0; i<ne; ++i)
		{
			FESurfaceElement& el = Element(i);
			if (m_bitfc && (el.m_elem[0] == nullptr)) FindElements(el);
			else if (el.m_elem[0] == nullptr) el.m_elem[0] = FindElement(el);
			//to make sure
			else if (m_bitfc && (el.m_elem[1] == nullptr)) FindElements(el);
			assert(el.m_elem[0] != nullptr);
		}
	}

	// serialize material point data
	for (int i = 0; i < Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		int nint = el.GaussPoints();
		for (int n = 0; n < nint; ++n)
		{
			FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(el.GetMaterialPoint(n));
			assert(pt);
			pt->Serialize(ar);
		}
	}
}

//-----------------------------------------------------------------------------
void FESurface::GetNodalCoordinates(FESurfaceElement& el, vec3d* rt)
{
	FEMesh& mesh = *GetMesh();
	int neln = el.Nodes();
    if (!m_bshellb) for (int j = 0; j < neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;
    else for (int j = 0; j < neln; ++j) rt[j] = mesh.Node(el.m_node[j]).st();
}

//-----------------------------------------------------------------------------
void FESurface::GetReferenceNodalCoordinates(FESurfaceElement& el, vec3d* r0)
{
	FEMesh& mesh = *GetMesh();
	int neln = el.Nodes();
    if (!m_bshellb) for (int j = 0; j < neln; ++j) r0[j] = mesh.Node(el.m_node[j]).m_r0;
    else for (int j = 0; j < neln; ++j) r0[j] = mesh.Node(el.m_node[j]).s0();
}

//-----------------------------------------------------------------------------
// Get current coordinates at intermediate configuration
void FESurface::GetNodalCoordinates(FESurfaceElement& el, double alpha, vec3d* rt)
{
	int neln = el.Nodes();
	for (int j = 0; j<neln; ++j) {
		FENode& node = Node(el.m_lnode[j]);
		rt[j] = node.m_rt*alpha + node.m_rp*(1.0 - alpha);
        if (m_bshellb) rt[j] -= node.m_dt*alpha + node.m_dp*(1.0 - alpha);
	}
}

//-----------------------------------------------------------------------------
// Evaluate field variables
double FESurface::Evaluate(FESurfaceMaterialPoint& mp, int dof)
{
	double v[FEElement::MAX_NODES];

	FESurfaceElement& el = *mp.SurfaceElement();
	int neln = el.Nodes();
	for (int j = 0; j < neln; ++j) v[j] = Node(el.m_lnode[j]).get(dof);

	double s = el.eval(v, mp.m_index);
	return s;
}

//-----------------------------------------------------------------------------
// Evaluate field variables
double FESurface::Evaluate(int nface, int dof)
{
    double v = 0;

    FESurfaceElement& el = Element(nface);
    int neln = el.Nodes();
    for (int j = 0; j < neln; ++j) v += Node(el.m_lnode[j]).get(dof);

    return v/neln;
}

//-----------------------------------------------------------------------------
void FESurface::LoadVector(FEGlobalVector& R, const FEDofList& dofList, bool breference, FESurfaceVectorIntegrand f)
{
	int dofPerNode = dofList.Size();
	int order = (dofPerNode == 1 ? dofList.InterpolationOrder(0) : -1);

	int NE = Elements();
	#pragma omp parallel for shared(R, dofList, f)
	for (int i = 0; i < NE; ++i)
	{
		vector<double> fe;
		vector<int> lm;
		vec3d re[FEElement::MAX_NODES];
		std::vector<double> G(dofPerNode, 0.0);

		// get the next element
		FESurfaceElement& el = Element(i);

		// init the element vector
		int neln = el.ShapeFunctions(order);
		int ndof = dofPerNode * neln;
		fe.assign(ndof, 0.0);

		// get the nodal coordinates
		if (breference)
			GetReferenceNodalCoordinates(el, re);
		else
			GetNodalCoordinates(el, re);

		// calculate element vector
		FESurfaceDofShape dof_a;
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n = 0; n < nint; ++n)
		{
			FESurfaceMaterialPoint& pt = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));

			// kinematics at integration points
			pt.dxr = el.eval_deriv1(re, n);
			pt.dxs = el.eval_deriv2(re, n);

			pt.m_shape = el.H(n);

			double* H = el.H(order, n);
			double* Hr = el.Gr(order, n);
			double* Hs = el.Gr(order, n);

			// put it all together
			for (int j = 0; j<neln; ++j)
			{
				// shape function and derivatives
				dof_a.index = j;
				dof_a.shape = H[j];
				dof_a.shape_deriv_r = Hr[j];
				dof_a.shape_deriv_s = Hs[j];

				// evaluate the integrand
				f(pt, dof_a, G);

				for (int k = 0; k < dofPerNode; ++k)
				{
					fe[dofPerNode * j + k] += G[k] * w[n];
				}
			}
		}

		// get the corresponding LM vector
		UnpackLM(el, dofList, lm);

		// Assemble into global vector
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FESurface::LoadStiffness(FELinearSystem& LS, const FEDofList& dofList_a, const FEDofList& dofList_b, FESurfaceMatrixIntegrand f)
{
	FEElementMatrix ke;

	int dofPerNode_a = dofList_a.Size();
	int dofPerNode_b = dofList_b.Size();

	int order_a = (dofPerNode_a == 1 ? dofList_a.InterpolationOrder(0) : -1);
	int order_b = (dofPerNode_b == 1 ? dofList_b.InterpolationOrder(0) : -1);

	vec3d rt[FEElement::MAX_NODES];

	matrix kab(dofPerNode_a, dofPerNode_b);
	FESurfaceDofShape dof_a, dof_b;

	int NE = Elements();
	for (int m = 0; m<NE; ++m)
	{
		// get the surface element
		FESurfaceElement& el = Element(m);

		ke.SetNodes(el.m_node);

		// shape functions
		int neln = el.Nodes();
		int nn_a = el.ShapeFunctions(dofPerNode_a);
		int nn_b = el.ShapeFunctions(dofPerNode_b);

		// get the element stiffness matrix
		int ndof_a = dofPerNode_a * nn_a;
		int ndof_b = dofPerNode_b * nn_b;
		ke.resize(ndof_a, ndof_b);

		// calculate element stiffness
		int nint = el.GaussPoints();

		// gauss weights
		double* w = el.GaussWeights();

		// nodal coordinates
		GetNodalCoordinates(el, rt);

		// repeat over integration points
		ke.zero();
		for (int n = 0; n<nint; ++n)
		{
			FESurfaceMaterialPoint& pt = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));

			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);

			// tangents at integration point
			pt.dxr = vec3d(0, 0, 0);
			pt.dxs = vec3d(0, 0, 0);
			for (int i = 0; i<neln; ++i)
			{
				pt.dxr += rt[i] * Gr[i];
				pt.dxs += rt[i] * Gs[i];
			}

			// calculate stiffness component
			for (int i = 0; i < nn_a; ++i)
			{
				// shape function values
				dof_a.index = i;
				dof_a.shape = el.H(order_a, n)[i];
				dof_a.shape_deriv_r = el.Gr(order_a, n)[i];
				dof_a.shape_deriv_s = el.Gs(order_a, n)[i];

				for (int j = 0; j < nn_b; ++j)
				{
					// shape function values
					dof_b.index = j;
					dof_b.shape = el.H(order_b, n)[j];
					dof_b.shape_deriv_r = el.Gr(order_b, n)[j];
					dof_b.shape_deriv_s = el.Gs(order_b, n)[j];

					// evaluate integrand
					kab.zero();
					f(pt, dof_a, dof_b, kab);

					// add it to the local element matrix
					ke.adds(dofPerNode_a * i, dofPerNode_b * j, kab, w[n]);
				}
			}
		}

		// get the element's LM vector
		std::vector<int>& lma = ke.RowIndices();
		std::vector<int>& lmb = ke.ColumnsIndices();
		UnpackLM(el, dofList_a, lma);
		UnpackLM(el, dofList_b, lmb);

		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}
