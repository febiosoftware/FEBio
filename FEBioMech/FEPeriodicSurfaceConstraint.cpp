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
#include "FEPeriodicSurfaceConstraint.h"
#include "FECore/FENormalProjection.h"
#include "FECore/FEGlobalMatrix.h"
#include <FECore/FELinearSystem.h>
#include "FECore/log.h"
#include <FECore/FEMesh.h>
#include "FEBioMech.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FEPeriodicSurfaceConstraint, FEContactInterface)
	ADD_PARAMETER(m_atol     , "tolerance");
	ADD_PARAMETER(m_eps      , "penalty");
	ADD_PARAMETER(m_btwo_pass, "two_pass");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Creates a surface for use with an FEPeriodicSurfaceConstraint interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FEPeriodicSurfaceConstraintSurface::Init()
{
	// always intialize base class first!
	if (FEContactSurface::Init() == false) return false;

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	m_gap.resize(nn);		// gap funtion
	m_pme.assign(nn, static_cast<FESurfaceElement*>(0));	// penetrated secondary element
	m_rs.resize(nn);		// natural coords of projected primary node on secondary element
	m_Lm.resize(nn);		// Lagrangian multipliers

	// set initial values
	zero(m_gap);
	zero(m_Lm);

	return true;
}

//-----------------------------------------------------------------------------
//! Calculate the center of mass for this surface
//!
vec3d FEPeriodicSurfaceConstraintSurface::CenterOfMass()
{
	vec3d c(0, 0, 0);
	int N = Nodes();
	for (int i = 0; i<N; ++i) c += Node(i).m_r0;
	if (N != 0) c /= (double)N;
	return c;
}

//-----------------------------------------------------------------------------
void FEPeriodicSurfaceConstraintSurface::Serialize(DumpStream& ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsShallow())
	{
		if (ar.IsSaving())
		{
			ar << m_Lm << m_gap;
		}
		else
		{
			ar >> m_Lm >> m_gap;
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			ar << m_gap;
			ar << m_rs;
			ar << m_Lm;
			ar << m_nref;
		}
		else
		{
			ar >> m_gap;
			ar >> m_rs;
			ar >> m_Lm;
			ar >> m_nref;
		}
	}
}

//-----------------------------------------------------------------------------
// FEPeriodicSurfaceConstraint
//-----------------------------------------------------------------------------

FEPeriodicSurfaceConstraint::FEPeriodicSurfaceConstraint(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem), m_dofU(pfem)
{
	static int count = 1;
	SetID(count++);

	m_stol = 0.01;
	m_srad = 1.0;
	m_atol = 0;
	m_eps = 0;
	m_btwo_pass = false;

	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}

	// set parents
	m_ss.SetContactInterface(this);
	m_ms.SetContactInterface(this);

	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------
bool FEPeriodicSurfaceConstraint::Init()
{
	// create the surfaces
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FEPeriodicSurfaceConstraint::BuildMatrixProfile(FEGlobalMatrix& K)
{
	FEMesh& mesh = GetMesh();

	// get the DOFS
	const int dof_X = GetDOFIndex("x");
	const int dof_Y = GetDOFIndex("y");
	const int dof_Z = GetDOFIndex("z");
	const int dof_RU = GetDOFIndex("Ru");
	const int dof_RV = GetDOFIndex("Rv");
	const int dof_RW = GetDOFIndex("Rw");

	vector<int> lm(6 * 5);

	int nref = m_ss.m_nref;
	FESurfaceElement* pref = m_ss.m_pme[nref];

	int n0 = pref->Nodes();
	int nr0[FEElement::MAX_NODES];
	for (int j = 0; j<n0; ++j) nr0[j] = pref->m_node[j];

	assign(lm, -1);

	lm[0] = m_ss.Node(nref).m_ID[dof_X];
	lm[1] = m_ss.Node(nref).m_ID[dof_Y];
	lm[2] = m_ss.Node(nref).m_ID[dof_Z];
	lm[3] = m_ss.Node(nref).m_ID[dof_RU];
	lm[4] = m_ss.Node(nref).m_ID[dof_RV];
	lm[5] = m_ss.Node(nref).m_ID[dof_RW];

	for (int k = 0; k<n0; ++k)
	{
		vector<int>& id = mesh.Node(nr0[k]).m_ID;
		lm[6 * (k + 1)] = id[dof_X];
		lm[6 * (k + 1) + 1] = id[dof_Y];
		lm[6 * (k + 1) + 2] = id[dof_Z];
		lm[6 * (k + 1) + 3] = id[dof_RU];
		lm[6 * (k + 1) + 4] = id[dof_RV];
		lm[6 * (k + 1) + 5] = id[dof_RW];
	}

	for (int j = 0; j<m_ss.Nodes(); ++j)
	{
		FESurfaceElement& me = *m_ss.m_pme[j];
		int* en = &me.m_node[0];

		assign(lm, -1);

		int n = me.Nodes();

		lm[0] = m_ss.Node(j).m_ID[dof_X];
		lm[1] = m_ss.Node(j).m_ID[dof_Y];
		lm[2] = m_ss.Node(j).m_ID[dof_Z];
		lm[3] = m_ss.Node(j).m_ID[dof_RU];
		lm[4] = m_ss.Node(j).m_ID[dof_RV];
		lm[5] = m_ss.Node(j).m_ID[dof_RW];

		for (int k = 0; k<n; ++k)
		{
			vector<int>& id = mesh.Node(en[k]).m_ID;
			lm[6 * (k + 1)] = id[dof_X];
			lm[6 * (k + 1) + 1] = id[dof_Y];
			lm[6 * (k + 1) + 2] = id[dof_Z];
			lm[6 * (k + 1) + 3] = id[dof_RU];
			lm[6 * (k + 1) + 4] = id[dof_RV];
			lm[6 * (k + 1) + 5] = id[dof_RW];
		}

		K.build_add(lm);
	}

}

//-----------------------------------------------------------------------------
void FEPeriodicSurfaceConstraint::Activate()
{
	// don't forget to call the base class
	FEContactInterface::Activate();

	// project primary surface onto secondary surface
	ProjectSurface(m_ss, m_ms, false);
	ProjectSurface(m_ms, m_ss, false);
}

//-----------------------------------------------------------------------------
//! project surface

void FEPeriodicSurfaceConstraint::ProjectSurface(FEPeriodicSurfaceConstraintSurface& ss, FEPeriodicSurfaceConstraintSurface& ms, bool bmove)
{
	FEMesh& mesh = GetMesh();

	FENormalProjection np(ms);
	np.SetTolerance(m_stol);
	np.SetSearchRadius(m_srad);
	np.Init();

	int i;
	double rs[2];

	// get the primary's center of mass
	vec3d cs = ss.CenterOfMass();

	// get the secondary's center of mass
	vec3d cm = ms.CenterOfMass();

	// get the relative distance
	vec3d cr = cm - cs;

	// unit vector in direction of cr
	// this will serve as the projection distance
	vec3d cn(cr); cn.unit();

	// loop over all primary nodes
	for (i = 0; i<ss.Nodes(); ++i)
	{
		FENode& node = ss.Node(i);

		// get the nodal position
		vec3d r0 = node.m_r0;

		// find the intersection with the secondary surface
		ss.m_pme[i] = np.Project(r0, cn, rs);
		assert(ss.m_pme[i]);

		ss.m_rs[i][0] = rs[0];
		ss.m_rs[i][1] = rs[1];
	}

	// if the reference node has not been located, we do it now.
	if (ss.m_nref < 0)
	{
		// we pick the node that is closest to the center of mass
		double dmin = (ss.Node(0).m_rt - cs).norm(), d;
		int nref = 0;
		for (i = 1; i<ss.Nodes(); ++i)
		{
			d = (ss.Node(i).m_rt - cs).norm();
			if (d < dmin)
			{
				dmin = d;
				nref = i;
			}
		}
		ss.m_nref = nref;
	}
}

//-----------------------------------------------------------------------------
void FEPeriodicSurfaceConstraint::Update()
{
	int i, j, ne, n0;
	FESurfaceElement* pme;

	FEMesh& mesh = *m_ss.GetMesh();

	vec3d us, um, u0;
	vec3d umi[FEElement::MAX_NODES];

	// update gap functions
	int npass = (m_btwo_pass ? 2 : 1);
	for (int np = 0; np<npass; ++np)
	{
		FEPeriodicSurfaceConstraintSurface& ss = (np == 0 ? m_ss : m_ms);
		FEPeriodicSurfaceConstraintSurface& ms = (np == 0 ? m_ms : m_ss);

		int N = ss.Nodes();

		// calculate the reference node displacement
		n0 = ss.m_nref;
		FENode& node = ss.Node(n0);
		us = node.m_rt - node.m_r0;

		// get the secondary element
		pme = ss.m_pme[n0];

		// calculate the secondary displacement
		ne = pme->Nodes();
		for (j = 0; j<ne; ++j)
		{
			FENode& node = ms.Node(pme->m_lnode[j]);
			umi[j] = node.m_rt - node.m_r0;
		}
		um = pme->eval(umi, ss.m_rs[n0][0], ss.m_rs[n0][1]);

		// calculate relative reference displacement
		u0 = us - um;

		for (i = 0; i<N; ++i)
		{
			// calculate the primary displacement
			FENode& node = ss.Node(i);
			us = node.m_rt - node.m_r0;

			// get the secondary element
			pme = ss.m_pme[i];

			// calculate the secondary displacement
			ne = pme->Nodes();
			for (j = 0; j<ne; ++j)
			{
				FENode& node = ms.Node(pme->m_lnode[j]);
				umi[j] = node.m_rt - node.m_r0;
			}
			um = pme->eval(umi, ss.m_rs[i][0], ss.m_rs[i][1]);

			// calculate gap function
			ss.m_gap[i] = us - um - u0;
		}
	}
}


//-----------------------------------------------------------------------------
void FEPeriodicSurfaceConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	int j, k, l, m, n;
	int nseln, nmeln;

	double *Gr, *Gs;

	// jacobian
	double detJ;

	vec3d dxr, dxs;
	vec3d rt[FEElement::MAX_NODES], r0[FEElement::MAX_NODES];
	double* w;

	// natural coordinates of primary node in secondary element
	double r, s;

	// contact force
	vec3d tc;

	// shape function values
	double N[FEElement::MAX_NODES];

	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	vector<int> sLM;
	vector<int> mLM;
	vector<int> LM0;

	int npass = (m_btwo_pass ? 2 : 1);
	for (int np = 0; np<npass; ++np)
	{
		FEPeriodicSurfaceConstraintSurface& ss = (np == 0 ? m_ss : m_ms);
		FEPeriodicSurfaceConstraintSurface& ms = (np == 0 ? m_ms : m_ss);

		// get the reference element
		int nref = ss.m_nref;
		FESurfaceElement* pref = ss.m_pme[nref];

		// loop over all primary facets
		int ne = ss.Elements();
		for (j = 0; j<ne; ++j)
		{
			// get the primary element
			FESurfaceElement& sel = ss.Element(j);

			// get the element's LM vector
			ss.UnpackLM(sel, sLM);

			nseln = sel.Nodes();

			for (int i = 0; i<nseln; ++i)
			{
				r0[i] = ss.GetMesh()->Node(sel.m_node[i]).m_r0;
				rt[i] = ss.GetMesh()->Node(sel.m_node[i]).m_rt;
			}

			w = sel.GaussWeights();

			// loop over primary element nodes (which are the integration points as well)
			for (n = 0; n<nseln; ++n)
			{
				Gr = sel.Gr(n);
				Gs = sel.Gs(n);

				m = sel.m_lnode[n];

				// calculate jacobian
				dxr = dxs = vec3d(0, 0, 0);
				for (k = 0; k<nseln; ++k)
				{
					dxr.x += Gr[k] * r0[k].x;
					dxr.y += Gr[k] * r0[k].y;
					dxr.z += Gr[k] * r0[k].z;

					dxs.x += Gs[k] * r0[k].x;
					dxs.y += Gs[k] * r0[k].y;
					dxs.z += Gs[k] * r0[k].z;
				}

				detJ = (dxr ^ dxs).norm();

				// get primary node contact force
				tc = ss.m_Lm[m] + ss.m_gap[m] * m_eps;

				// get the secondary element
				FESurfaceElement& mel = *ss.m_pme[m];
				ms.UnpackLM(mel, mLM);

				nmeln = mel.Nodes();

				// isoparametric coordinates of the projected primary node
				// onto the secondary element
				r = ss.m_rs[m][0];
				s = ss.m_rs[m][1];

				// get the secondary shape function values at this primary node
				if (nmeln == 4)
				{
					// quadrilateral
					N[0] = 0.25*(1 - r)*(1 - s);
					N[1] = 0.25*(1 + r)*(1 - s);
					N[2] = 0.25*(1 + r)*(1 + s);
					N[3] = 0.25*(1 - r)*(1 + s);
				}
				else if (nmeln == 3)
				{
					// triangle
					N[0] = 1 - r - s;
					N[1] = r;
					N[2] = s;
				}

				// calculate force vector
				fe.resize(3 * (nmeln + 1));
				fe[0] = -detJ*w[n] * tc.x;
				fe[1] = -detJ*w[n] * tc.y;
				fe[2] = -detJ*w[n] * tc.z;
				for (l = 0; l<nmeln; ++l)
				{
					fe[3 * (l + 1)] = detJ*w[n] * tc.x*N[l];
					fe[3 * (l + 1) + 1] = detJ*w[n] * tc.y*N[l];
					fe[3 * (l + 1) + 2] = detJ*w[n] * tc.z*N[l];
				}

				// fill the lm array
				lm.resize(3 * (nmeln + 1));
				lm[0] = sLM[n * 3];
				lm[1] = sLM[n * 3 + 1];
				lm[2] = sLM[n * 3 + 2];
				for (l = 0; l<nmeln; ++l)
				{
					lm[3 * (l + 1)] = mLM[l * 3];
					lm[3 * (l + 1) + 1] = mLM[l * 3 + 1];
					lm[3 * (l + 1) + 2] = mLM[l * 3 + 2];
				}

				// fill the en array
				en.resize(nmeln + 1);
				en[0] = sel.m_node[n];
				for (l = 0; l<nmeln; ++l) en[l + 1] = mel.m_node[l];

				// assemble into global force vector
				R.Assemble(en, lm, fe);
			}
		}
	}
}


//-----------------------------------------------------------------------------
void FEPeriodicSurfaceConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	int j, k, l, n, m;
	int nseln, nmeln, ndof;

	FEElementMatrix ke;

	vector<int> lm(15);
	vector<int> en(5);

	double *Gr, *Gs, *w;

	vec3d rt[FEElement::MAX_NODES], r0[FEElement::MAX_NODES];
	vec3d rtm[FEElement::MAX_NODES];

	double detJ, r, s;
	vec3d dxr, dxs;
	double H[FEElement::MAX_NODES];

	vec3d gap, Lm, tc;

	// curvature tensor K
	double K[2][2] = { 0 };

	//	double scale = -0.0035*m_fem.GetMesh().GetBoundingBox().radius();

	vector<int> sLM;
	vector<int> mLM;
	vector<int> LM0;

	int n0[FEElement::MAX_NODES];

	int npass = (m_btwo_pass ? 2 : 1);
	for (int np = 0; np<npass; ++np)
	{
		FEPeriodicSurfaceConstraintSurface& ss = (np == 0 ? m_ss : m_ms);
		FEPeriodicSurfaceConstraintSurface& ms = (np == 0 ? m_ms : m_ss);

		// get the reference element
		int nref = ss.m_nref;
		FESurfaceElement* pref = ss.m_pme[nref];

		// grab the data we'll need for this element
		ms.UnpackLM(*pref, LM0);
		int ne0 = pref->Nodes();
		for (j = 0; j<ne0; ++j) n0[j] = pref->m_node[j];
		r = ss.m_rs[nref][0];
		s = ss.m_rs[nref][1];
		double N0[FEElement::MAX_NODES];
		pref->shape_fnc(N0, r, s);

		// number of degrees of freedom
		ndof = 3 * (1 + ne0);
		vector<double> N(ndof / 3);
		N[0] = 1;
		for (k = 0; k<ne0; ++k) N[k + 1] = -N0[k];

		// stiffness matrix
		ke.resize(ndof, ndof); ke.zero();
		for (k = 0; k<ndof / 3; ++k)
			for (l = 0; l<ndof / 3; ++l)
			{
				ke[3 * k][3 * l] = m_eps*N[k] * N[l];
				ke[3 * k + 1][3 * l + 1] = m_eps*N[k] * N[l];
				ke[3 * k + 2][3 * l + 2] = m_eps*N[k] * N[l];
			}

		// fill the lm array
		lm.resize(3 * (ne0 + 1));
		lm[0] = ss.Node(nref).m_ID[m_dofU[0]];
		lm[1] = ss.Node(nref).m_ID[m_dofU[1]];
		lm[2] = ss.Node(nref).m_ID[m_dofU[2]];
		for (l = 0; l<ne0; ++l)
		{
			lm[3 * (l + 1)] = LM0[l * 3];
			lm[3 * (l + 1) + 1] = LM0[l * 3 + 1];
			lm[3 * (l + 1) + 2] = LM0[l * 3 + 2];
		}

		// fill the en array
		en.resize(ne0 + 1);
		en[0] = ss.NodeIndex(nref);
		for (l = 0; l<ne0; ++l) en[l + 1] = n0[l];

		// assemble stiffness matrix
		//		psolver->AssembleStiffness(en, lm, ke);

		// loop over all primary elements
		int ne = ss.Elements();
		for (j = 0; j<ne; ++j)
		{
			FESurfaceElement& se = ss.Element(j);

			// get the element's LM vector
			ss.UnpackLM(se, sLM);

			nseln = se.Nodes();

			for (int i = 0; i<nseln; ++i)
			{
				r0[i] = ss.GetMesh()->Node(se.m_node[i]).m_r0;
				rt[i] = ss.GetMesh()->Node(se.m_node[i]).m_rt;
			}

			w = se.GaussWeights();

			// loop over all integration points (that is nodes)
			for (n = 0; n<nseln; ++n)
			{
				Gr = se.Gr(n);
				Gs = se.Gs(n);

				m = se.m_lnode[n];

				// calculate jacobian
				dxr = dxs = vec3d(0, 0, 0);
				for (k = 0; k<nseln; ++k)
				{
					dxr.x += Gr[k] * r0[k].x;
					dxr.y += Gr[k] * r0[k].y;
					dxr.z += Gr[k] * r0[k].z;

					dxs.x += Gs[k] * r0[k].x;
					dxs.y += Gs[k] * r0[k].y;
					dxs.z += Gs[k] * r0[k].z;
				}

				detJ = (dxr ^ dxs).norm();

				// get the secondary element
				FESurfaceElement& me = *ss.m_pme[m];
				ms.UnpackLM(me, mLM);

				nmeln = me.Nodes();

				// get the secondary element node positions
				for (k = 0; k<nmeln; ++k) rtm[k] = ms.GetMesh()->Node(me.m_node[k]).m_rt;

				// primary node natural coordinates in secondary element
				r = ss.m_rs[m][0];
				s = ss.m_rs[m][1];

				// get primary node normal force
				tc = ss.m_Lm[m] + ss.m_gap[m] * m_eps; //ss.T[m];

				// get the secondary shape function values at this primary node
				if (nmeln == 4)
				{
					// quadrilateral
					H[0] = 0.25*(1 - r)*(1 - s);
					H[1] = 0.25*(1 + r)*(1 - s);
					H[2] = 0.25*(1 + r)*(1 + s);
					H[3] = 0.25*(1 - r)*(1 + s);
				}
				else if (nmeln == 3)
				{
					// triangle
					H[0] = 1 - r - s;
					H[1] = r;
					H[2] = s;
				}

				// number of degrees of freedom
				ndof = 3 * (1 + nmeln);

				vector<double> N(ndof / 3);
				N[0] = 1;
				for (k = 0; k<nmeln; ++k) N[k + 1] = -H[k];

				ke.resize(ndof, ndof); ke.zero();
				for (k = 0; k<ndof / 3; ++k)
					for (l = 0; l<ndof / 3; ++l)
					{
						ke[3 * k][3 * l] = w[n] * detJ*m_eps*N[k] * N[l];
						ke[3 * k + 1][3 * l + 1] = w[n] * detJ*m_eps*N[k] * N[l];
						ke[3 * k + 2][3 * l + 2] = w[n] * detJ*m_eps*N[k] * N[l];
					}

				// fill the lm array
				lm.resize(3 * (nmeln + 1));
				lm[0] = sLM[n * 3];
				lm[1] = sLM[n * 3 + 1];
				lm[2] = sLM[n * 3 + 2];
				for (l = 0; l<nmeln; ++l)
				{
					lm[3 * (l + 1)] = mLM[l * 3];
					lm[3 * (l + 1) + 1] = mLM[l * 3 + 1];
					lm[3 * (l + 1) + 2] = mLM[l * 3 + 2];
				}

				// fill the en array
				en.resize(nmeln + 1);
				en[0] = se.m_node[n];
				for (l = 0; l<nmeln; ++l) en[l + 1] = me.m_node[l];

				// assemble stiffness matrix
				ke.SetNodes(en);
				ke.SetIndices(lm);
				LS.Assemble(ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEPeriodicSurfaceConstraint::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != 1) return true;

	int i;
	bool bconv = true;

	double g;
	vec3d lm;

	// calculate initial norms
	double normL0 = 0;
	for (i = 0; i<m_ss.Nodes(); ++i)
	{
		lm = m_ss.m_Lm[i];
		normL0 += lm*lm;
	}
	for (i = 0; i<m_ms.Nodes(); ++i)
	{
		lm = m_ms.m_Lm[i];
		normL0 += lm*lm;
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i = 0; i<m_ss.Nodes(); ++i)
	{
		lm = m_ss.m_Lm[i] + m_ss.m_gap[i] * m_eps;

		normL1 += lm*lm;
		g = m_ss.m_gap[i].norm();
		normgc += g*g;
		++N;
	}

	for (i = 0; i<m_ms.Nodes(); ++i)
	{
		lm = m_ms.m_Lm[i] + m_ms.m_gap[i] * m_eps;

		normL1 += lm*lm;
		g = m_ms.m_gap[i].norm();
		normgc += g*g;
		++N;
	}
	if (N == 0) N = 1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// check convergence of constraints
	feLog(" surface constraint# %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0) / normL1);
	feLog("    normal force : %15le %15le\n", pctn, m_atol);
	feLog("    gap function : %15le       ***\n", normgc);

	if (pctn >= m_atol)
	{
		bconv = false;
		for (i = 0; i<m_ss.Nodes(); ++i)
		{
			// update Lagrange multipliers
			m_ss.m_Lm[i] = m_ss.m_Lm[i] + m_ss.m_gap[i] * m_eps;
		}

		for (i = 0; i<m_ms.Nodes(); ++i)
		{
			// update Lagrange multipliers
			m_ms.m_Lm[i] = m_ms.m_Lm[i] + m_ms.m_gap[i] * m_eps;
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
void FEPeriodicSurfaceConstraint::Serialize(DumpStream &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);
}
