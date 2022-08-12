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
#include "FERigidWallInterface.h"
#include <FECore/FENNQuery.h>
#include <FECore/FEModel.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/log.h>
#include <FEBioMech/FEElasticShellDomainOld.h>
#include <FECore/FELinearSystem.h>

// Macauley bracket
#define MBRACKET(x) ((x)>=0? (x): 0)
#define HEAVYSIDE(x) ((x)>=0?1:0)

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FERigidWallInterface, FESurfaceConstraint)
	ADD_PARAMETER(m_laugon , "laugon"    )->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
	ADD_PARAMETER(m_atol   , "tolerance" );
	ADD_PARAMETER(m_eps    , "penalty"   );
	ADD_PARAMETER(m_d	   , "offset"    );
	ADD_PARAMETER(m_a      , 4, "plane"  );
END_FECORE_CLASS();

///////////////////////////////////////////////////////////////////////////////
// FERigidWallSurface
///////////////////////////////////////////////////////////////////////////////


FERigidWallSurface::FERigidWallSurface(FEModel* pfem) : FESurface(pfem) 
{ 
	m_NQ.Attach(this); 

	// I want to use the FEModel class for this, but don't know how
	DOFS& dofs = pfem->GetDOFS();
	m_dofX = dofs.GetDOF("x");
	m_dofY = dofs.GetDOF("y");
	m_dofZ = dofs.GetDOF("z");
}

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FERigidWallSurface::Init()
{
	// always intialize base class first!
	if (FESurface::Init() == false) return false;

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	m_gap.assign(nn, 0);		// gap funtion
	m_nu.resize(nn);		// node normal 
	m_pme.assign(nn, static_cast<FESurfaceElement*>(0));		// penetrated secondary surface element
	m_rs.resize(nn);		// natural coords of projected primary node on secondary element
	m_rsp.resize(nn);
	m_Lm.assign(nn, 0);
	m_M.resize(nn);
	m_Lt.resize(nn);
	m_off.assign(nn, 0.0);
	m_eps.assign(nn, 1.0);

	// we calculate the gap offset values
	// This value is used to take the shell thickness into account
	// note that we force rigid shells to have zero thickness
	FEMesh& m = *m_pMesh;
	vector<double> tag(m.Nodes());
	zero(tag);
	for (int nd=0; nd<m.Domains(); ++nd)
	{
		FEElasticShellDomainOld* psd = dynamic_cast<FEElasticShellDomainOld*>(&m.Domain(nd));
		if (psd)
		{
			for (int i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				int n = el.Nodes();
				for (int j=0; j<n; ++j) tag[el.m_node[j]] = 0.5*el.m_h0[j];
			}
		}
	}
	for (int i=0; i<nn; ++i) m_off[i] = tag[NodeIndex(i)];

	return true;
}


//-----------------------------------------------------------------------------
//! 

vec3d FERigidWallSurface::traction(int inode)
{
	vec3d t(0,0,0);
	FESurfaceElement* pe = m_pme[inode];
	if (pe)
	{
		FESurfaceElement& el = *pe;
		double Tn = m_Lm[inode];
		double T1 = m_Lt[inode][0];
		double T2 = m_Lt[inode][1];
		double r = m_rs[inode][0];
		double s = m_rs[inode][1];

		vec3d tn = m_nu[inode]*Tn, tt;
		vec3d e[2];
		ContraBaseVectors0(el, r, s, e);
		tt = e[0]*T1 + e[1]*T2;
		t = tn + tt;
	}

	return t;
}

//-----------------------------------------------------------------------------

void FERigidWallSurface::UpdateNormals()
{
	int i, j, jp1, jm1;
	int N = Nodes();
	int NE = Elements();
	for (i=0; i<N; ++i) m_nu[i] = vec3d(0,0,0);
	vec3d y[FEElement::MAX_NODES], e1, e2;

	for (i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		int ne = el.Nodes();
		for (j=0; j<ne; ++j) y[j] = Node(el.m_lnode[j]).m_rt;

		for (j=0; j<ne; ++j)
		{
			jp1 = (j+1)%ne;
			jm1 = (j+ne-1)%ne;

			e1 = y[jp1] - y[j];
			e2 = y[jm1] - y[j];

			m_nu[el.m_lnode[j]] -= e1 ^ e2;						
		}
	}

	for (i=0; i<N; ++i) m_nu[i].unit();
}

//-----------------------------------------------------------------------------
void FERigidWallSurface::Serialize(DumpStream &ar)
{
	FESurface::Serialize(ar);
	ar & m_gap;
	ar & m_nu;
	ar & m_rs;
	ar & m_rsp;
	ar & m_Lm;
	ar & m_M;
	ar & m_Lt;
	ar & m_off;
	ar & m_eps;
}

//-----------------------------------------------------------------------------
void FERigidWallSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];
	}
}

///////////////////////////////////////////////////////////////////////////////
// FERigidWallInterface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FERigidWallInterface::FERigidWallInterface(FEModel* pfem) : FESurfaceConstraint(pfem)
{
	static int count = 1;
	SetID(count++);

	m_laugon = 0;

	m_a[0] = m_a[1] = m_a[2] = m_a[3] = 0.0;

	m_ss = new FERigidWallSurface(pfem);
	m_eps = 0;
	m_atol = 0;
	m_d = 0.0;
};

//-----------------------------------------------------------------------------
FERigidWallInterface::~FERigidWallInterface()
{
	delete m_ss;
}

//-----------------------------------------------------------------------------
//! Initializes the rigid wall interface data

bool FERigidWallInterface::Init()
{
	// create the surface
	if (m_ss->Init() == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FERigidWallInterface::BuildMatrixProfile(FEGlobalMatrix& K)
{
	FERigidWallSurface& surf = *m_ss;

	FEModel& fem = *GetFEModel();

	// get the DOFS
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");
	const int dof_RU = fem.GetDOFIndex("Ru");
	const int dof_RV = fem.GetDOFIndex("Rv");
	const int dof_RW = fem.GetDOFIndex("Rw");

	vector<int> lm(6);
	for (int j=0; j< surf.Nodes(); ++j)
	{
		if (surf.m_gap[j] >= 0)
		{
			lm[0] = surf.Node(j).m_ID[dof_X];
			lm[1] = surf.Node(j).m_ID[dof_Y];
			lm[2] = surf.Node(j).m_ID[dof_Z];
			lm[3] = surf.Node(j).m_ID[dof_RU];
			lm[4] = surf.Node(j).m_ID[dof_RV];
			lm[5] = surf.Node(j).m_ID[dof_RW];

			K.build_add(lm);
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidWallInterface::Activate()
{
	// don't forget to call the base class
	FESurfaceConstraint::Activate();

	// project primary surface onto secondary surface
	ProjectSurface(*m_ss);
}

//-----------------------------------------------------------------------------
//!  Projects the primary surface onto the plane

void FERigidWallInterface::ProjectSurface(FERigidWallSurface& ss)
{
	FERigidWallSurface& surf = *m_ss;

	// loop over all primary surface nodes
	for (int i=0; i< surf.Nodes(); ++i)
	{
		// get the nodal position
		vec3d r = surf.Node(i).m_rt;

		// project this node onto the plane
		vec3d q = ProjectToPlane(r);

		// get the local surface normal
		vec3d np = PlaneNormal(q);

		// calculate offset
		q += np*m_d;

		// the normal is set to the secondary surface element normal
		surf.m_nu[i] = np;
	
		// calculate initial gap
		surf.m_gap[i] = -(np*(r - q)) + surf.m_off[i];
	}
}

//-----------------------------------------------------------------------------
//!  Updates rigid wall data

void FERigidWallInterface::Update()
{
	// project primary surface onto secondary surface
	ProjectSurface(*m_ss);
}

//-----------------------------------------------------------------------------

void FERigidWallInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	int j, k, m, n;
	int nseln;

	double *Gr, *Gs;

	// jacobian
	double detJ;

	vec3d dxr, dxs;
	vec3d rt[FEElement::MAX_NODES], r0[FEElement::MAX_NODES];
	double* w;

	// normal force
	double tn;

	// element contact force vector
	// note that this assumes that the element to which the integration node
	// connects is a four noded quadrilateral
	vector<double> fe(3);

	// the lm array for this force vector
	vector<int> lm(3);

	// the en array
	vector<int> en(1);

	vector<int> sLM;

	// penalty value
	double pen = m_eps, eps;
	
	// loop over all primary surface facets
	FERigidWallSurface& surf = *m_ss;

	int ne = surf.Elements();
	for (j=0; j<ne; ++j)
	{
		// get the next element
		FESurfaceElement& sel = surf.Element(j);

		// get the element's LM vector
		surf.UnpackLM(sel, sLM);

		nseln = sel.Nodes();

		for (int i=0; i<nseln; ++i)
		{
			r0[i] = surf.GetMesh()->Node(sel.m_node[i]).m_r0;
			rt[i] = surf.GetMesh()->Node(sel.m_node[i]).m_rt;
		}
		w = sel.GaussWeights();

		// loop over element nodes (which are the integration points as well)
		for (n=0; n<nseln; ++n)
		{
			Gr = sel.Gr(n);
			Gs = sel.Gs(n);

			m = sel.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a secondary surface element associated with it
//			if (ss.Lm[m] >= 0)
			{
				// calculate jacobian
				dxr = dxs = vec3d(0,0,0);
				for (k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				detJ = (dxr ^ dxs).norm();

				// get node normal force
				eps = pen* surf.m_eps[m];
				tn = surf.m_Lm[m] + eps* surf.m_gap[m];
				tn = MBRACKET(tn);

				// get the node normal
				vec3d& nu = surf.m_nu[m];

				// calculate force vector
				fe[0] = detJ*w[n]*tn*nu.x;
				fe[1] = detJ*w[n]*tn*nu.y;
				fe[2] = detJ*w[n]*tn*nu.z;
	
				// fill the lm array
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				// fill the en array
				en[0] = sel.m_node[n];

				// assemble into global force vector
				R.Assemble(en, lm, fe);
			}
		}
	}

}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness contribution for the rigid wall
//! interface.
//! \todo I think there are a couple of stiffness terms missing in this formulation

void FERigidWallInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	int j, k, l, n, m;
	int nseln, ndof;

	FEElementMatrix ke;

	vector<int> lm(3);
	vector<int> en(1);

	double *Gr, *Gs, *w;
	vec3d rt[FEElement::MAX_NODES], r0[FEElement::MAX_NODES];

	double detJ, tn;
	vec3d dxr, dxs;
	double N[3];

	double gap, Lm;

	vector<int> sLM;

	// penalty value
	double pen = m_eps, eps;

	// loop over all primary surface elements
	FERigidWallSurface& surf = *m_ss;
	int ne = surf.Elements();
	for (j=0; j<ne; ++j)
	{
		FESurfaceElement& se = surf.Element(j);

		// get the element's LM vector
		surf.UnpackLM(se, sLM);

		nseln = se.Nodes();

		for (int i=0; i<nseln; ++i)
		{
			r0[i] = surf.GetMesh()->Node(se.m_node[i]).m_r0;
			rt[i] = surf.GetMesh()->Node(se.m_node[i]).m_rt;
		}

		w = se.GaussWeights();

		// loop over all integration points (that is nodes)
		for (n=0; n<nseln; ++n)
		{
			Gr = se.Gr(n);
			Gs = se.Gs(n);

			m = se.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a secondary surface element associated with it
//			if (ss.Lm[m] >= 0)
			{
				// calculate jacobian
				dxr = dxs = vec3d(0,0,0);
				for (k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				detJ = (dxr ^ dxs).norm();

				// gap
				gap = surf.m_gap[m];

				// lagrange multiplier
				Lm = surf.m_Lm[m];

				// get node normal force
				eps = pen* surf.m_eps[m];

				tn = surf.m_Lm[m] + eps* surf.m_gap[m];
				tn = MBRACKET(tn);

				// get the node normal
				vec3d& nu = surf.m_nu[m];

				// set up the N vector
				N[0] = nu.x;
				N[1] = nu.y;
				N[2] = nu.z;

				ndof = 3;

				// fill stiffness matrix
				// TODO: I don't think this is correct, since
				// if the rigid wall is not a plance, some terms
				// are probably missing
				ke.resize(ndof, ndof);
				for (k=0; k<ndof; ++k)
					for (l=0; l<ndof; ++l)
						ke[k][l] = w[n]*detJ*eps*HEAVYSIDE(Lm+eps*gap)*N[k]*N[l];

				// create lm array
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				// create the en array
				en[0] = se.m_node[n];
						
				// assemble stiffness matrix
				ke.SetNodes(en);
				ke.SetIndices(lm);
				LS.Assemble(ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------

bool FERigidWallInterface::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != 1) return true;

	int i;
	double Lm;
	bool bconv = true;

	FERigidWallSurface& surf = *m_ss;

	// penalty value
	double pen = m_eps, eps;

	// calculate initial norms
	double normL0 = 0;
	for (i=0; i< surf.Nodes(); ++i) normL0 += surf.m_Lm[i]* surf.m_Lm[i];
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i=0; i< surf.Nodes(); ++i)
	{
		// update Lagrange multipliers
		eps = pen* surf.m_eps[i];

		Lm = surf.m_Lm[i] + eps* surf.m_gap[i];
		Lm = MBRACKET(Lm);
		normL1 += Lm*Lm;
		if (surf.m_gap[i] > 0)
		{
			normgc += surf.m_gap[i]* surf.m_gap[i];
			++N;
		}
	}	
	if (N==0) N = 1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// check convergence of constraints
	feLog(" rigid wall interface # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	feLog("    normal force : %15le %15le\n", pctn, m_atol);
	feLog("    gap function : %15le       ***\n", normgc);
		
	if (pctn >= m_atol)
	{
		bconv = false;
		for (i=0; i< surf.Nodes(); ++i)
		{
			// update Lagrange multipliers
			eps = pen* surf.m_eps[i];

			Lm = surf.m_Lm[i] + eps* surf.m_gap[i];
			surf.m_Lm[i] = MBRACKET(Lm);
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------
vec3d FERigidWallInterface::PlaneNormal(const vec3d& r)
{
	vec3d n(m_a[0], m_a[1], m_a[2]);
	n.unit();
	return n;
}

//-----------------------------------------------------------------------------
vec3d FERigidWallInterface::ProjectToPlane(const vec3d& r)
{
	double d = m_a[3];

	double l = m_a[0] * r.x + m_a[1] * r.y + m_a[2] * r.z - d;
	return vec3d(r.x - l * m_a[0], r.y - l * m_a[1], r.z - l * m_a[2]);
}

//-----------------------------------------------------------------------------
void FERigidWallInterface::Serialize(DumpStream &ar)
{
	FESurfaceConstraint::Serialize(ar);
	m_ss->Serialize(ar);
}
