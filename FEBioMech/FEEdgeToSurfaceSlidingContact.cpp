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
#include "FEEdgeToSurfaceSlidingContact.h"
#include <FECore/FENode.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FEBox.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEClosestPointProjection.h>
#include <FECore/FEGlobalMatrix.h>
#include <stdexcept>

void FEEdgeToSurfaceSlidingContactSurface::Update()
{
	for (int i = 0; i < Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);

		vec3d re[FEElement::MAX_NODES];
		GetNodalCoordinates(el, re);

		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FESurfaceMaterialPoint& mp = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));
			mp.m_rt = el.eval(re, n);
			
			// kinematics at integration points
			mp.dxr = el.eval_deriv1(re, n);
			mp.dxs = el.eval_deriv2(re, n);

			mp.m_Jt = (mp.dxr ^ mp.dxs).norm();
		}
	}
}

FEEdgeToSurfaceSlidingContactSurface::FEEdgeToSurfaceSlidingContactSurface(FEModel* fem) : FEContactSurface(fem)
{

}

//=================================================================================================
FEEdgeToSurfaceSlidingContactEdge::FEEdgeToSurfaceSlidingContactEdge(FEModel* fem) : FEEdge(fem)
{

}

bool FEEdgeToSurfaceSlidingContactEdge::Init()
{
	if (!FEEdge::Init()) return false;
	DOFS& dofs = GetFEModel()->GetDOFS();
	m_dofX = dofs.GetDOF("x");
	m_dofY = dofs.GetDOF("y");
	m_dofZ = dofs.GetDOF("z");
	m_points.resize(Nodes());
	return true;
}

void FEEdgeToSurfaceSlidingContactEdge::UnpackLM(FELineElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N * 3);
	for (int i = 0; i < N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		lm[3 * i    ] = id[m_dofX];
		lm[3 * i + 1] = id[m_dofY];
		lm[3 * i + 2] = id[m_dofZ];
	}
}

FEMaterialPoint* FEEdgeToSurfaceSlidingContactEdge::CreateMaterialPoint()
{
	return new FEE2SSlidingContactPoint();
}

bool FEEdgeToSurfaceSlidingContactEdge::Create(FESegmentSet& eset)
{
	return FEEdge::Create(eset, FE_LINE2NI);
}

void FEEdgeToSurfaceSlidingContactEdge::Update()
{
	for (int i = 0; i < Elements(); ++i)
	{
		FELineElement& el = Element(i);

		vec3d rt[FEElement::MAX_NODES];
		GetNodalCoordinates(el, rt);

		// Jacobian is length in spatial coordinates divided by length in parametric coordinates (=2)
		double J = (rt[1] - rt[0]).Length() / 2.0;

		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FEE2SSlidingContactPoint& mp = static_cast<FEE2SSlidingContactPoint&>(*el.GetMaterialPoint(n));
			mp.m_rt = el.eval(rt, n);
			mp.m_Jt = J;
		}
	}
}

void FEEdgeToSurfaceSlidingContactEdge::Serialize(DumpStream& ar)
{
	FEEdge::Serialize(ar);
	ar.LockPointerTable();
	for (int i = 0; i < m_points.size(); ++i)
	{
		ar & m_points[i];
	}
	ar.UnlockPointerTable();
}

//=================================================================================================
BEGIN_FECORE_CLASS(FEEdgeToSurfaceSlidingContact, FESurfaceConstraint)
	ADD_PARAMETER(m_atol         , "tolerance"    );
	ADD_PARAMETER(m_eps          , "penalty"      );
	ADD_PARAMETER(m_naugmin      , "minaug"       );
	ADD_PARAMETER(m_naugmax      , "maxaug"       );
	ADD_PARAMETER(m_stol         , "search_tol"   );
	ADD_PARAMETER(m_nsegup       , "seg_up"       );
	ADD_PARAMETER(m_sradius      , "search_radius");

	ADD_PROPERTY(m_edge, "edgelist")->AddFlag(FEProperty::Reference);

END_FECORE_CLASS();

FEEdgeToSurfaceSlidingContact::FEEdgeToSurfaceSlidingContact(FEModel* fem) : FESurfaceConstraint(fem), m_edge(fem), m_surf(fem)
{
	m_atol = 0.1;
	m_eps = 1.0;
	m_naugmin = 0;
	m_naugmax = 0;
	m_stol = 0.01;
	m_nsegup = 0;
	m_sradius = 0;

	m_bfirst = true;
}

FESurface* FEEdgeToSurfaceSlidingContact::GetSurface()
{
	return &m_surf;
}

// initialization
bool FEEdgeToSurfaceSlidingContact::Init()
{
	if (FESurfaceConstraint::Init() == false) return false;

	m_bfirst = true;

	if (m_edge.Init() == false) return false;
	if (m_surf.Init() == false) return false;

	return true;
}

//!  Projects the primary surface onto the secondary surface.
//!  That is, for each primary surface node we determine the closest
//!  secondary surface element and the projection of that node onto
//!  this element.

//! \todo this function needs to identify the different types of contact:
//!   1/ first contact
//!   2/ crossing of element boundary
//!	  3/ contact termination 
//!			either by failure to find projection or when g < tolerance

void FEEdgeToSurfaceSlidingContact::ProjectSurface(bool bupseg, bool bmove)
{
	FEClosestPointProjection cpp(m_surf);
	cpp.SetTolerance(m_stol);
	cpp.SetSearchRadius(m_sradius);
	cpp.HandleSpecialCases(true);
	cpp.Init();

	for (int i = 0; i < m_edge.Nodes(); ++i)
	{
		// get the node
		FENode& node = m_edge.Node(i);

		// get the nodal position
		vec3d x = node.m_rt;

		// get the global node number
		int m = m_edge.NodeIndex(i);

		FEE2SSlidingContactPoint& cp = m_edge.m_points[i];

		// get the previous secondary surface element (if any)
		FESurfaceElement* pme = cp.m_pme;

		// If the node is in contact, let's see if the node still is 
		// on the same element
		vec3d q;
		if (pme != 0)
		{
			FESurfaceElement& mel = *pme;

			double r = cp.m_rs[0];
			double s = cp.m_rs[1];

			q = m_surf.ProjectToSurface(mel, x, r, s);
			cp.m_rs[1] = s;
			cp.m_rs[0] = r;

			// we only check when we can update the segments
			// otherwise, we just stick with this element, even
			// if the node is no longer inside it.
			if (bupseg)
			{
				if (!m_surf.IsInsideElement(mel, r, s, m_stol))
				{
					// see if the node might have moved to another element
					FESurfaceElement* pold = pme;
					cp.m_rs = vec2d(0, 0);

					pme = cpp.Project(m, q, cp.m_rs);
				}
			}
		}
		else if (bupseg)
		{
			// get the secondary surface element
			// don't forget to initialize the search for the first node!
			cp.m_rs = vec2d(0, 0);
			pme = cpp.Project(m, q, cp.m_rs);
		}

		// if we found a secondary surface element, update the gap and normal data
		cp.m_pme = pme;
		if (pme != 0)
		{
			FESurfaceElement& mel = *cp.m_pme;

			double r = cp.m_rs[0];
			double s = cp.m_rs[1];

			// the normal is set to the secondary surface element normal
			cp.m_nu = m_surf.SurfaceNormal(mel, r, s);

			// calculate gap
			cp.m_gap = -(cp.m_nu * (x - q));
			if (bmove && (cp.m_gap > 0))
			{
				node.m_r0 = node.m_rt = q;
				cp.m_gap = 0;
			}

			// TODO: what should we do if the gap function becomes
			// negative? setting the Lagrange multipliers to zero
			// might make the system unstable.
/*			if (ss.gap[i] < 0)
			{
				ss.Lm[i] = 0;
				ss.pme[i] = 0;
			}
*/
		}
		else
		{
			// TODO: Is this a good criteria for out-of-contact?
			//		 perhaps this is not even necessary.
			// since the node is not in contact, we set the gap function 
			// and Lagrangian multiplier to zero
			cp.m_gap = 0;
			cp.m_Lm = 0;
		}
	}
}

// update
void FEEdgeToSurfaceSlidingContact::Update()
{
	FESurfaceConstraint::Update();

	int niter = GetFEModel()->GetCurrentStep()->GetFESolver()->m_niter;

	// should we do a segment update or not?
	// TODO: check what happens when m_nsegup == -1 and m_npass = 2;
	// We have to make sure that in this case, both surfaces get at least
	// one pass!
	bool bupdate = (m_bfirst || (m_nsegup == 0) ? true : (niter <= m_nsegup));

	m_edge.Update();
	m_surf.Update();

	// project primary surface onto secondary surface
	// this also calculates the nodal gap functions
	ProjectSurface(bupdate);

	// set the first-entry-flag to false
	m_bfirst = false;
}

// Build the matrix profile
void FEEdgeToSurfaceSlidingContact::BuildMatrixProfile(FEGlobalMatrix& K)
{
	// TODO: this is currently for max 6 nodes (hence 7=6+1)
	vector<int> lm(6 * 7);

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the DOFS
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");
	const int dof_RU = fem.GetDOFIndex("Ru");
	const int dof_RV = fem.GetDOFIndex("Rv");
	const int dof_RW = fem.GetDOFIndex("Rw");

	for (int j = 0; j < m_edge.Nodes(); ++j)
	{
		FESurfaceElement* pe = m_edge.m_points[j].m_pme;
		if (pe != 0)
		{
			FESurfaceElement& me = *pe;
			int* en = &me.m_node[0];

			// Note that we need to grab the rigid degrees of freedom as well
			// this is in case one of the nodes belongs to a rigid body.
			int n = me.Nodes();
			if (n == 3)
			{
				lm[6 * (3 + 1)    ] = -1; lm[6 * (3 + 2)    ] = -1; lm[6 * (3 + 3)    ] = -1;
				lm[6 * (3 + 1) + 1] = -1; lm[6 * (3 + 2) + 1] = -1; lm[6 * (3 + 3) + 1] = -1;
				lm[6 * (3 + 1) + 2] = -1; lm[6 * (3 + 2) + 2] = -1; lm[6 * (3 + 3) + 2] = -1;
				lm[6 * (3 + 1) + 3] = -1; lm[6 * (3 + 2) + 3] = -1; lm[6 * (3 + 3) + 3] = -1;
				lm[6 * (3 + 1) + 4] = -1; lm[6 * (3 + 2) + 4] = -1; lm[6 * (3 + 3) + 4] = -1;
				lm[6 * (3 + 1) + 5] = -1; lm[6 * (3 + 2) + 5] = -1; lm[6 * (3 + 3) + 5] = -1;
			}
			if (n == 4)
			{
				lm[6 * (4 + 1)    ] = -1; lm[6 * (4 + 2)    ] = -1;
				lm[6 * (4 + 1) + 1] = -1; lm[6 * (4 + 2) + 1] = -1;
				lm[6 * (4 + 1) + 2] = -1; lm[6 * (4 + 2) + 2] = -1;
				lm[6 * (4 + 1) + 3] = -1; lm[6 * (4 + 2) + 3] = -1;
				lm[6 * (4 + 1) + 4] = -1; lm[6 * (4 + 2) + 4] = -1;
				lm[6 * (4 + 1) + 5] = -1; lm[6 * (4 + 2) + 5] = -1;
			}

			lm[0] = m_edge.Node(j).m_ID[dof_X];
			lm[1] = m_edge.Node(j).m_ID[dof_Y];
			lm[2] = m_edge.Node(j).m_ID[dof_Z];
			lm[3] = m_edge.Node(j).m_ID[dof_RU];
			lm[4] = m_edge.Node(j).m_ID[dof_RV];
			lm[5] = m_edge.Node(j).m_ID[dof_RW];

			for (int k = 0; k < n; ++k)
			{
				vector<int>& id = mesh.Node(en[k]).m_ID;
				lm[6 * (k + 1)    ] = id[dof_X];
				lm[6 * (k + 1) + 1] = id[dof_Y];
				lm[6 * (k + 1) + 2] = id[dof_Z];
				lm[6 * (k + 1) + 3] = id[dof_RU];
				lm[6 * (k + 1) + 4] = id[dof_RV];
				lm[6 * (k + 1) + 5] = id[dof_RW];
			}

			K.build_add(lm);
		}
	}
}

// The LoadVector function evaluates the "forces" that contribute to the residual of the system
void FEEdgeToSurfaceSlidingContact::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	// the elements LM vectors
	vector<int> sLM;
	vector<int> mLM;

	const int MN = FEElement::MAX_NODES;
	vec3d r0[MN];
	double w[MN];
	double* Gr, * Gs;
	double detJ[MN];
	vec3d dxr, dxs;

	int ne = m_edge.Elements();
	for (int j = 0; j < ne; ++j)
	{
		// get the next element
		FELineElement& sel = m_edge.Element(j);
		int nseln = sel.Nodes();

		// get the element's LM array
		m_edge.UnpackLM(sel, sLM);

		// nodal coordinates
		for (int i = 0; i < nseln; ++i) r0[i] = m_edge.GetMesh()->Node(sel.m_node[i]).m_r0;

		// we calculate all the metrics we need before we
		// calculate the nodal forces
		double* w1 = sel.GaussWeights();
		for (int n = 0; n < nseln; ++n)
		{
			FEE2SSlidingContactPoint& mp1 = static_cast<FEE2SSlidingContactPoint&>(*sel.GetMaterialPoint(n));
			detJ[n] = mp1.m_Jt;
			w[n] = sel.GaussWeights()[n];
		}

		// loop over primary surface element nodes (which are the integration points as well)
		// and calculate the contact nodal force
		for (int n = 0; n < nseln; ++n)
		{
			// get the local node number
			int m = sel.m_lnode[n];

			FEE2SSlidingContactPoint& cp = m_edge.m_points[m];

			// see if this node's constraint is active
			// that is, if it has an element associated with it
			// TODO: is this a good way to test for an active constraint
			// The rigid wall criteria seems to work much better.
			if (cp.m_pme != 0)
			{
				// This node is active and could lead to a non-zero
				// contact force.
				// get the secondary surface element
				FESurfaceElement& mel = *cp.m_pme;
				m_surf.UnpackLM(mel, mLM);

				// calculate the degrees of freedom
				int nmeln = mel.Nodes();
				int ndof = 3 * (nmeln + 1);
				fe.resize(ndof);

				// calculate the nodal force
				ContactNodalForce(m, mel, fe);

				// multiply force with weights
				for (int l = 0; l < ndof; ++l) fe[l] *= detJ[n] * w[n];

				// fill the lm array
				lm.resize(3 * (nmeln + 1));
				lm[0] = sLM[n * 3    ];
				lm[1] = sLM[n * 3 + 1];
				lm[2] = sLM[n * 3 + 2];

				for (int l = 0; l < nmeln; ++l)
				{
					lm[3 * (l + 1)    ] = mLM[l * 3];
					lm[3 * (l + 1) + 1] = mLM[l * 3 + 1];
					lm[3 * (l + 1) + 2] = mLM[l * 3 + 2];
				}

				// fill the en array
				en.resize(nmeln + 1);
				en[0] = sel.m_node[n];
				for (int l = 0; l < nmeln; ++l) en[l + 1] = mel.m_node[l];
					
				// assemble into global force vector
				R.Assemble(en, lm, fe);
			}
		}
	}
}

void FEEdgeToSurfaceSlidingContact::ContactNodalForce(int m, FESurfaceElement& mel, vector<double>& fe)
{
	// max nr of element nodes
	const int MAXMN = FEElement::MAX_NODES;

	// secondary surface element nodes
	vec3d rtm[MAXMN];

	// shape function values
	double H[MAXMN], Hr[MAXMN], Hs[MAXMN];

	// contact vectors
	double N[3 * (MAXMN + 1)], N1[3 * (MAXMN + 1)], N2[3 * (MAXMN + 1)];
	double T1[3 * (MAXMN + 1)], T2[3 * (MAXMN + 1)], D1[3 * (MAXMN + 1)], D2[3 * (MAXMN + 1)];

	// surface metrics
	double A[2][2], M[2][2], K[2][2];
	double detA;

	double eps = m_eps;

	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	double Tt[2];

	int nmeln, ndof;

	FEE2SSlidingContactPoint& cp = m_edge.m_points[m];

	// gap function
	double gap = cp.m_gap;

	// get primary surface node normal force
	double Ln = cp.m_Lm;
	double tn = Ln + eps * gap;
	tn = MBRACKET(tn);

	// get the primary surface  node normal
	vec3d& nu = cp.m_nu;

	nmeln = mel.Nodes();
	ndof = 3 * (1 + nmeln);

	// get the secondary surface element node positions
	for (int k = 0; k < nmeln; ++k) rtm[k] = mesh.Node(mel.m_node[k]).m_rt;

	// isoparametric coordinates of the projected node
	// onto the secondary surface element
	double r = cp.m_rs[0];
	double s = cp.m_rs[1];

	// get the secondary surface element shape function values at this node
	mel.shape_fnc(H, r, s);

	// --- N O R M A L   T R A C T I O N ---

	// calculate contact vectors for normal traction
	N[0] = nu.x;
	N[1] = nu.y;
	N[2] = nu.z;
	for (int l = 0; l < nmeln; ++l)
	{
		N[3 * (l + 1)] = -H[l] * nu.x;
		N[3 * (l + 1) + 1] = -H[l] * nu.y;
		N[3 * (l + 1) + 2] = -H[l] * nu.z;
	}

	// calculate force vector
	for (int l = 0; l < ndof; ++l) fe[l] = tn * N[l];
}

// Evaluates the contriubtion to the stiffness matrix
void FEEdgeToSurfaceSlidingContact::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEElementMatrix ke;

	const int MAXMN = FEElement::MAX_NODES;
	vector<int> lm(3 * (MAXMN + 1));
	vector<int> en(MAXMN + 1);

	double* Gr, * Gs, w[MAXMN];
	vec3d r0[MAXMN];

	double detJ[MAXMN];
	vec3d dxr, dxs;

	vector<int> sLM;
	vector<int> mLM;

	// loop over all primary surface elements
	int ne = m_edge.Elements();
	for (int j = 0; j < ne; ++j)
	{
		// unpack the next element
		FELineElement& se = m_edge.Element(j);
		int nseln = se.Nodes();

		// get the element's LM array
		m_edge.UnpackLM(se, sLM);

		// get the nodal coordinates
		for (int i = 0; i < nseln; ++i) r0[i] = m_edge.GetMesh()->Node(se.m_node[i]).m_r0;

		// get all the metrics we need 
		for (int n = 0; n < nseln; ++n)
		{
			FEE2SSlidingContactPoint& mp1 = static_cast<FEE2SSlidingContactPoint&>(*se.GetMaterialPoint(n));
			detJ[n] = mp1.m_Jt;
			w[n] = se.GaussWeights()[n];
		}

		// loop over all integration points (that is nodes)
		for (int n = 0; n < nseln; ++n)
		{
			int m = se.m_lnode[n];

			FEE2SSlidingContactPoint& cp = m_edge.m_points[m];

			// see if this node's constraint is active
			// that is, if it has an element associated with it
			if (cp.m_pme != 0)
			{
				// get the secondary surface element
				FESurfaceElement& me = *cp.m_pme;

				// get the secondary surface element's LM array
				m_surf.UnpackLM(me, mLM);

				int nmeln = me.Nodes();
				int ndof = 3 * (nmeln + 1);

				// calculate the stiffness matrix
				ke.resize(ndof, ndof);
				ContactNodalStiffness(m, me, ke);

				// muliply with weights
				for (int k = 0; k < ndof; ++k)
					for (int l = 0; l < ndof; ++l) ke[k][l] *= detJ[n] * w[n];

				// fill the lm array
				lm[0] = sLM[n * 3    ];
				lm[1] = sLM[n * 3 + 1];
				lm[2] = sLM[n * 3 + 2];

				for (int k = 0; k < nmeln; ++k)
				{
					lm[3 * (k + 1)    ] = mLM[k * 3];
					lm[3 * (k + 1) + 1] = mLM[k * 3 + 1];
					lm[3 * (k + 1) + 2] = mLM[k * 3 + 2];
				}

				// create the en array
				en.resize(nmeln + 1);
				en[0] = se.m_node[n];
				for (int k = 0; k < nmeln; ++k) en[k + 1] = me.m_node[k];

				// assemble stiffness matrix
				ke.SetNodes(en);
				ke.SetIndices(lm);
				LS.Assemble(ke);
			}
		}
	}
}

void FEEdgeToSurfaceSlidingContact::ContactNodalStiffness(int m, FESurfaceElement& mel, matrix& ke)
{
	const int MAXMN = FEElement::MAX_NODES;

	vector<int> lm(3 * (MAXMN + 1));
	vector<int> en(MAXMN + 1);

	vec3d dxr, dxs;
	double H[MAXMN], Hr[MAXMN], Hs[MAXMN];

	double N[3 * (MAXMN + 1)], T1[3 * (MAXMN + 1)], T2[3 * (MAXMN + 1)];
	double N1[3 * (MAXMN + 1)], N2[3 * (MAXMN + 1)], D1[3 * (MAXMN + 1)], D2[3 * (MAXMN + 1)];
	double Nb1[3 * (MAXMN + 1)], Nb2[3 * (MAXMN + 1)];

	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// nr of element nodes and degrees of freedom 
	int nmeln = mel.Nodes();
	int ndof = 3 * (1 + nmeln);

	// penalty factor
	double eps = m_eps;

	// nodal coordinates
	vec3d rt[MAXMN];
	for (int j = 0; j < nmeln; ++j) rt[j] = mesh.Node(mel.m_node[j]).m_rt;

	// node natural coordinates in secondary surface element
	FEE2SSlidingContactPoint& cp = m_edge.m_points[m];
	double r = cp.m_rs[0];
	double s = cp.m_rs[1];

	// gap
	double gap = cp.m_gap;

	// lagrange multiplier
	double Lm = cp.m_Lm;

	// get node normal force
	double tn = Lm + eps * gap;
	tn = MBRACKET(tn);

	// get the node normal
	vec3d& nu = cp.m_nu;

	// get the secondary surface element shape function values and the derivatives at this node
	mel.shape_fnc(H, r, s);
	mel.shape_deriv(Hr, Hs, r, s);

	// get the tangent vectors
	vec3d tau[2];
	m_surf.CoBaseVectors(mel, r, s, tau);

	// set up the N vector
	N[0] = nu.x;
	N[1] = nu.y;
	N[2] = nu.z;

	for (int k = 0; k < nmeln; ++k)
	{
		N[(k + 1) * 3] = -H[k] * nu.x;
		N[(k + 1) * 3 + 1] = -H[k] * nu.y;
		N[(k + 1) * 3 + 2] = -H[k] * nu.z;
	}

	// set up the Ti vectors
	T1[0] = tau[0].x; T2[0] = tau[1].x;
	T1[1] = tau[0].y; T2[1] = tau[1].y;
	T1[2] = tau[0].z; T2[2] = tau[1].z;

	for (int k = 0; k < nmeln; ++k)
	{
		T1[(k + 1) * 3    ] = -H[k] * tau[0].x;
		T1[(k + 1) * 3 + 1] = -H[k] * tau[0].y;
		T1[(k + 1) * 3 + 2] = -H[k] * tau[0].z;

		T2[(k + 1) * 3    ] = -H[k] * tau[1].x;
		T2[(k + 1) * 3 + 1] = -H[k] * tau[1].y;
		T2[(k + 1) * 3 + 2] = -H[k] * tau[1].z;
	}

	// set up the Ni vectors
	N1[0] = N2[0] = 0;
	N1[1] = N2[1] = 0;
	N1[2] = N2[2] = 0;

	for (int k = 0; k < nmeln; ++k)
	{
		N1[(k + 1) * 3    ] = -Hr[k] * nu.x;
		N1[(k + 1) * 3 + 1] = -Hr[k] * nu.y;
		N1[(k + 1) * 3 + 2] = -Hr[k] * nu.z;

		N2[(k + 1) * 3    ] = -Hs[k] * nu.x;
		N2[(k + 1) * 3 + 1] = -Hs[k] * nu.y;
		N2[(k + 1) * 3 + 2] = -Hs[k] * nu.z;
	}

	// calculate metric tensor
	mat2d M;
	M[0][0] = tau[0] * tau[0]; M[0][1] = tau[0] * tau[1];
	M[1][0] = tau[1] * tau[0]; M[1][1] = tau[1] * tau[1];

	// calculate reciprocal metric tensor
	mat2d Mi = M.inverse();

	// calculate curvature tensor
	double K[2][2] = { 0 };
	double Grr[FEElement::MAX_NODES];
	double Grs[FEElement::MAX_NODES];
	double Gss[FEElement::MAX_NODES];
	mel.shape_deriv2(Grr, Grs, Gss, r, s);
	for (int k = 0; k < nmeln; ++k)
	{
		K[0][0] += (nu * rt[k]) * Grr[k];
		K[0][1] += (nu * rt[k]) * Grs[k];
		K[1][0] += (nu * rt[k]) * Grs[k];
		K[1][1] += (nu * rt[k]) * Gss[k];
	}

	// setup A matrix A = M + gK
	double A[2][2];
	A[0][0] = M[0][0] + gap * K[0][0];
	A[0][1] = M[0][1] + gap * K[0][1];
	A[1][0] = M[1][0] + gap * K[1][0];
	A[1][1] = M[1][1] + gap * K[1][1];

	// calculate determinant of A
	double detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	// setup Di vectors
	for (int k = 0; k < ndof; ++k)
	{
		D1[k] = (1 / detA) * (A[1][1] * (T1[k] + gap * N1[k]) - A[0][1] * (T2[k] + gap * N2[k]));
		D2[k] = (1 / detA) * (A[0][0] * (T2[k] + gap * N2[k]) - A[0][1] * (T1[k] + gap * N1[k]));
	}

	// setup Nbi vectors
	for (int k = 0; k < ndof; ++k)
	{
		Nb1[k] = N1[k] - K[0][1] * D2[k];
		Nb2[k] = N2[k] - K[0][1] * D1[k];
	}

	// --- N O R M A L   S T I F F N E S S ---
	double sum;
	for (int k = 0; k < ndof; ++k)
		for (int l = 0; l < ndof; ++l)
		{
			sum = 0;

			sum = Mi[0][0] * Nb1[k] * Nb1[l] + Mi[0][1] * (Nb1[k] * Nb2[l] + Nb2[k] * Nb1[l]) + Mi[1][1] * Nb2[k] * Nb2[l];
			sum *= gap;
			sum -= D1[k] * N1[l] + D2[k] * N2[l] + N1[k] * D1[l] + N2[k] * D2[l];
			sum += K[0][1] * (D1[k] * D2[l] + D2[k] * D1[l]);
			sum *= tn;

			sum += eps * HEAVYSIDE(Lm + eps * gap) * N[k] * N[l];

			ke[k][l] = sum;
		}
}

void FEEdgeToSurfaceSlidingContact::Serialize(DumpStream& ar)
{
	FESurfaceConstraint::Serialize(ar);
	m_edge.Serialize(ar);
	m_surf.Serialize(ar);
}
