/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FERigidSlidingContact.h"
#include <FECore/FEModel.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/log.h>
#include "FEMechModel.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FERigidSlidingContact, FEContactInterface)
	ADD_PARAMETER(m_atol     , "tolerance");
	ADD_PARAMETER(m_eps      , "penalty"  );
	ADD_PARAMETER(m_gtol     , "gaptol"   );
	ADD_PARAMETER(m_naugmin  , "minaug"   );
	ADD_PARAMETER(m_naugmax  , "maxaug"   );
	ADD_PARAMETER(m_bautopen , "auto_penalty");
	ADD_PARAMETER(m_rigidName, "rigid");
END_FECORE_CLASS();

///////////////////////////////////////////////////////////////////////////////
// FERigidSphereSurface
///////////////////////////////////////////////////////////////////////////////

void FERigidSlidingSurface::DATA::Serialize(DumpStream& ar)
{
	ar & gap;
	ar & nu;
	ar & Lm;
	ar & eps;
}


FERigidSlidingSurface::FERigidSlidingSurface(FEModel* pfem) : FEContactSurface(pfem)
{ 
	m_NQ.Attach(this); 

	m_Fc = vec3d(0,0,0);

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

bool FERigidSlidingSurface::Init()
{
	// always intialize base class first!
	if (FESurface::Init() == false) return false;

	// get the number of elements
	int NE = Elements();

	// count total number of integration points
	int nintTotal = 0;
	for (int i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		nintTotal += el.GaussPoints();
	}

	// allocate other surface data
	m_data.resize(nintTotal);
	for (int i=0; i<nintTotal; ++i)
	{
		DATA& d = m_data[i];
		d.gap = 0.0;
		d.Lm = 0.0;
		d.nu = vec3d(0,0,0);
		d.eps = 1.0;
	}

	return true;
}


//-----------------------------------------------------------------------------
// TODO: I don't think we need this
vec3d FERigidSlidingSurface::traction(int inode)
{
	return vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
vec3d FERigidSlidingSurface::GetContactForce()
{
	return m_Fc;
}

//-----------------------------------------------------------------------------
void FERigidSlidingSurface::Serialize(DumpStream &ar)
{
	FESurface::Serialize(ar);
	ar & m_data & m_Fc;
}

//-----------------------------------------------------------------------------
void FERigidSlidingSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_lnode[i];
		FENode& node = Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];
	}
}

///////////////////////////////////////////////////////////////////////////////
// FERigidSlidingContact
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FERigidSlidingContact::FERigidSlidingContact(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem)
{
	static int count = 1;
	SetID(count++);

	m_rigid = 0;
	m_eps = 0;
	m_atol = 0;
	m_gtol = 0;
	m_naugmin = 0;
	m_naugmax = 10;
	m_bautopen = false;
	m_rigidName[0] = 0;
};

//-----------------------------------------------------------------------------
FERigidSlidingContact::~FERigidSlidingContact()
{
	m_rigid = 0;
}

//-----------------------------------------------------------------------------
//! Initializes the rigid wall interface data

bool FERigidSlidingContact::Init()
{
	// make sure a rigid surface was defined
	if (m_rigid == 0)
	{
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		m_rigid = fem.FindRigidSurface(m_rigidName);
		if (m_rigid == 0) return false;
	}

	// create the surface
	if (m_ss.Init() == false) return false;

	// initialize rigid surface
	if (m_rigid->Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FERigidSlidingContact::CalcAutoPenalty(FERigidSlidingSurface& s)
{
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();

	// loop over all surface elements
	int c = 0;
	for (int i = 0; i<s.Elements(); ++i)
	{
		// get the surface element
		FESurfaceElement& el = s.Element(i);

		// calculate a modulus
		double eps = AutoPenalty(el, s);

		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j = 0; j<nint; ++j)
		{
			FERigidSlidingSurface::DATA& pt = s.m_data[c++];
			pt.eps = eps;
		}
	}
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FERigidSlidingContact::BuildMatrixProfile(FEGlobalMatrix& K)
{
	FEModel& fem = *GetFEModel();

	// get the DOFS
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");
	const int dof_RU = fem.GetDOFIndex("Ru");
	const int dof_RV = fem.GetDOFIndex("Rv");
	const int dof_RW = fem.GetDOFIndex("Rw");

	int c = 0;
	vector<int> lm;
	lm.reserve(3*FEElement::MAX_NODES);
	for (int i=0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& el = m_ss.Element(i);

		// see if any integration points are in contact
		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FERigidSlidingSurface::DATA& d = m_ss.m_data[c + j];
			if (d.gap >= 0)
			{
				m_ss.UnpackLM(el, lm);
				K.build_add(lm);
				break;
			}
		}
		c += el.GaussPoints();
	}
}

//-----------------------------------------------------------------------------
void FERigidSlidingContact::Activate()
{
	// don't forget to call the base class
	FEContactInterface::Activate();

	// calculate penalty factors
	if (m_bautopen) CalcAutoPenalty(m_ss);

	// project primary surface onto secondary surface
	ProjectSurface(m_ss);
}

//-----------------------------------------------------------------------------
//!  Projects the primary surface onto the plane

void FERigidSlidingContact::ProjectSurface(FERigidSlidingSurface& ss)
{
	vec3d rt[FEElement::MAX_NODES];

	// loop over all primary surface elements
	int c = 0;
	for (int i=0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& el = m_ss.Element(i);
		int neln = el.Nodes();
		int nint = el.GaussPoints();

		// get the nodal coordinates
		for (int j=0; j<neln; ++j) rt[j] = m_ss.Node(el.m_lnode[j]).m_rt;

		// loop over all integration points
		for (int j=0; j<nint; ++j)
		{
			// get the integration point data
			FERigidSlidingSurface::DATA& d = m_ss.m_data[c++];

			// get the nodal position
			vec3d r = el.Evaluate(rt, j);

			// project this node onto the plane
			vec3d q = m_rigid->Project(r);

			// get the local surface normal
			vec3d np = m_rigid->Normal(q);

			// the normal is set to the secondary surface element normal
			d.nu = np;
	
			// calculate initial gap
			d.gap = -(np*(r - q));
		}
	}
}

//-----------------------------------------------------------------------------
//!  Updates rigid wall data

void FERigidSlidingContact::Update()
{
	// project primary surface onto secondary surface
	ProjectSurface(m_ss);
}

//-----------------------------------------------------------------------------

void FERigidSlidingContact::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	vector<int> lm;
	const int MELN = FEElement::MAX_NODES;
	double detJ[MELN], w[MELN];
	vec3d r0[MELN];
	vector<double> fe;

	// zero total force
	m_ss.m_Fc = vec3d(0,0,0);

	// loop over all primary surface elements
	int c = 0;
	for (int i = 0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& se = m_ss.Element(i);
		int neln = se.Nodes();
		int nint = se.GaussPoints();
		int ndof = 3*neln;

		// get the element's LM vector
		m_ss.UnpackLM(se, lm);

		// nodal coordinates
		for (int j = 0; j<neln; ++j) r0[j] = m_ss.Node(se.m_lnode[j]).m_r0;

		// we calculate all the metrics we need before we
		// calculate the nodal forces
		for (int j = 0; j<nint; ++j)
		{
			double* Gr = se.Gr(j);
			double* Gs = se.Gs(j);

			// calculate jacobian
			// note that we are integrating over the reference surface
			vec3d dxr, dxs;
			for (int k = 0; k<neln; ++k)
			{
				dxr.x += Gr[k] * r0[k].x;
				dxr.y += Gr[k] * r0[k].y;
				dxr.z += Gr[k] * r0[k].z;

				dxs.x += Gs[k] * r0[k].x;
				dxs.y += Gs[k] * r0[k].y;
				dxs.z += Gs[k] * r0[k].z;
			}

			// jacobians
			detJ[j] = (dxr ^ dxs).norm();

			// integration weights
			w[j] = se.GaussWeights()[j];
		}

		// loop over all integration points
		for (int j = 0; j<nint; ++j)
		{
			// get integration point data
			FERigidSlidingSurface::DATA& d = m_ss.m_data[c++];

			// calculate shape functions
			double* H = se.H(j);

			// get normal vector
			vec3d nu = d.nu;

			// gap function
			double g = d.gap;

			// lagrange multiplier
			double Lm = d.Lm;

			// penalty value
			double eps = m_eps*d.eps;

			// contact traction
			double tn = Lm + eps*g;
			tn = MBRACKET(tn);

			// calculate the force vector
			fe.resize(ndof);

			for (int k = 0; k<neln; ++k)
			{
				fe[3 * k    ] = H[k] * nu.x;
				fe[3 * k + 1] = H[k] * nu.y;
				fe[3 * k + 2] = H[k] * nu.z;
			}
			for (int k = 0; k<ndof; ++k) fe[k] *= tn*detJ[j] * w[j];

			// add to total reaction force
			for (int k=0; k<neln; ++k)
			{
				m_ss.m_Fc.x += fe[3*k  ];
				m_ss.m_Fc.y += fe[3*k+1];
				m_ss.m_Fc.z += fe[3*k+2];
			}

			// assemble the global residual
			R.Assemble(se.m_node, lm, fe);
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidSlidingContact::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	vector<int> lm;
	const int MELN = FEElement::MAX_NODES;
	const int MINT = FEElement::MAX_INTPOINTS;
	vec3d r0[MELN];
	double detJ[MINT], w[MINT];
	double N[3*MELN];
	FEElementMatrix ke;

	// loop over all primary surface elements
	int c = 0;
	for (int i = 0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& se = m_ss.Element(i);
		int neln = se.Nodes();
		int nint = se.GaussPoints();
		int ndof = 3*neln;

		// get the element's LM vector
		m_ss.UnpackLM(se, lm);

		// nodal coordinates
		for (int j = 0; j<neln; ++j) r0[j] = m_ss.Node(se.m_lnode[j]).m_r0;

		// we calculate all the metrics we need before we
		// calculate the nodal forces
		for (int j = 0; j<nint; ++j)
		{
			double* Gr = se.Gr(j);
			double* Gs = se.Gs(j);

			// calculate jacobian
			// note that we are integrating over the reference surface
			vec3d dxr, dxs;
			for (int k = 0; k<neln; ++k)
			{
				dxr.x += Gr[k] * r0[k].x;
				dxr.y += Gr[k] * r0[k].y;
				dxr.z += Gr[k] * r0[k].z;

				dxs.x += Gs[k] * r0[k].x;
				dxs.y += Gs[k] * r0[k].y;
				dxs.z += Gs[k] * r0[k].z;
			}

			// jacobians
			detJ[j] = (dxr ^ dxs).norm();

			// integration weights
			w[j] = se.GaussWeights()[j];
		}

		// loop over all integration points
		for (int j = 0; j<nint; ++j)
		{
			// get integration point data
			FERigidSlidingSurface::DATA& d = m_ss.m_data[c++];

			// calculate shape functions
			double* H = se.H(j);

			// get normal vector
			vec3d nu = d.nu;

			// gap function
			double g = d.gap;

			// when the node is on the surface, the gap value
			// can flip-flop between positive and negative.
			if (fabs(g)<1e-20) g = 0;

			// lagrange multiplier
			double Lm = d.Lm;

			// penalty value
			double eps = m_eps*d.eps;

			// contact traction
			double tn = Lm + eps*g;
			tn = MBRACKET(tn);

			double dtn = eps*HEAVYSIDE(Lm + eps*g);

			// calculate the N-vector
			for (int k = 0; k<neln; ++k)
			{
				N[3 * k    ] = H[k] * nu.x;
				N[3 * k + 1] = H[k] * nu.y;
				N[3 * k + 2] = H[k] * nu.z;
			}

			// create the stiffness matrix
			ke.resize(ndof, ndof);

			// add the first order term (= D(tn)*dg )
			for (int k = 0; k<ndof; ++k)
				for (int l = 0; l<ndof; ++l) ke[k][l] = dtn*N[k] * N[l] * detJ[j] * w[j];

			// assemble the global residual
			ke.SetNodes(se.m_node);
			ke.SetIndices(lm);
			LS.Assemble(ke);
		}
	}
}

//-----------------------------------------------------------------------------
bool FERigidSlidingContact::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != 1) return true;

	// calculate initial norms
	double normL0 = 0;
	int c = 0;
	for (int i=0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& el = m_ss.Element(i);
		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FERigidSlidingSurface::DATA& d = m_ss.m_data[c++];
			normL0 += d.Lm*d.Lm;
		}
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	c = 0;
	for (int i = 0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& el = m_ss.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FERigidSlidingSurface::DATA& d = m_ss.m_data[c++];
			
			// update Lagrange multipliers
			double eps = m_eps*d.eps;
			double Lm = d.Lm + eps*d.gap;
			Lm = MBRACKET(Lm);
			normL1 += Lm*Lm;
			if (d.gap > 0)
			{
				normgc += d.gap*d.gap;
				++N;
			}
		}	
	}
	if (N==0) N = 1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// calculate and print convergence norms
	double lnorm = 0;
	if (normL1 != 0) lnorm = fabs(normL1 - normL0) / normL1; else lnorm = fabs(normL1 - normL0);

	// check convergence of constraints
	feLog(" rigid sliding contact # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	feLog("    normal force : %15le", lnorm);
	if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
	feLog("    gap function : %15le", normgc);
	if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");

	// check convergence
	bool bconv = true;
	if ((m_atol > 0) && (lnorm >= m_atol)) bconv = false;
	if ((m_gtol > 0) && (normgc > m_gtol)) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if (m_naugmax <= naug) bconv = true;

	if (bconv == false)
	{
		bconv = false;
		c = 0;
		for (int i = 0; i<m_ss.Elements(); ++i)
		{
			FESurfaceElement& el = m_ss.Element(i);
			for (int j = 0; j<el.GaussPoints(); ++j)
			{
				FERigidSlidingSurface::DATA& d = m_ss.m_data[c++];
	
				double eps = m_eps*d.eps;
				double Lm = d.Lm + eps*d.gap;
				d.Lm = MBRACKET(Lm);
			}	
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------

void FERigidSlidingContact::Serialize(DumpStream &ar)
{
	FEContactInterface::Serialize(ar);
	m_ss.Serialize(ar);
	if (m_rigid) m_rigid->Serialize(ar);
}
