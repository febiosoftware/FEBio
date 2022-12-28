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
#include "FETiedInterface.h"
#include <FECore/FEClosestPointProjection.h>
#include <FECore/FELinearSystem.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedInterface, FEContactInterface)
	ADD_PARAMETER(m_atol    , "tolerance"       );
	ADD_PARAMETER(m_eps     , "penalty"         );
	ADD_PARAMETER(m_naugmin , "minaug"          );
	ADD_PARAMETER(m_naugmax , "maxaug"          );
	ADD_PARAMETER(m_stol    , "search_tolerance");
	ADD_PARAMETER(m_boffset , "offset_shells"   );
	ADD_PARAMETER(m_Dmax    , "max_distance"    );
	ADD_PARAMETER(m_bspecial, "special"         );
	ADD_PARAMETER(m_breloc  , "node_reloc"      );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. Initialize default values.
FETiedInterface::FETiedInterface(FEModel* pfem) : FEContactInterface(pfem), ss(pfem), ms(pfem)
{
	static int count = 1;
	SetID(count++);

	// define sibling relationships
	ss.SetSibling(&ms);
	ms.SetSibling(&ss);

	// initial parameter values
	m_laugon = 0;
	m_atol = 0.01;
	m_eps = 1.0;
	m_stol = 0.0001;
	m_naugmin = 0;
	m_naugmax = 10;
	m_boffset = false;
	m_Dmax = 0.0;
	m_bspecial = true;
	m_breloc = false;

	// set parents
	ss.SetContactInterface(this);
	ms.SetContactInterface(this);
}

//-----------------------------------------------------------------------------
//! Initialization. This function intializes the surfaces data and projects the
//! primary surface onto the secondary surface.
//! 
bool FETiedInterface::Init()
{
	// set surface options
	ss.SetShellOffset(m_boffset);

	// create the surfaces
	if (ss.Init() == false) return false;
	if (ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FETiedInterface::BuildMatrixProfile(FEGlobalMatrix& K)
{
	FEMesh& mesh = GetMesh();

	// get the DOFS
	const int dof_X = GetDOFIndex("x");
	const int dof_Y = GetDOFIndex("y");
	const int dof_Z = GetDOFIndex("z");
	const int dof_RU = GetDOFIndex("Ru");
	const int dof_RV = GetDOFIndex("Rv");
	const int dof_RW = GetDOFIndex("Rw");

	if (m_laugon != 2)
	{
		const int LMSIZE = 6 * (FEElement::MAX_NODES + 1);
		vector<int> lm(LMSIZE);

		for (int j = 0; j < ss.Nodes(); ++j)
		{
			FESurfaceElement* pe = ss.m_data[j].m_pme;
			if (pe != 0)
			{
				FESurfaceElement& me = *pe;
				int* en = &me.m_lnode[0];

				int n = me.Nodes();
				lm.assign(LMSIZE, -1);

				lm[0] = ss.Node(j).m_ID[dof_X];
				lm[1] = ss.Node(j).m_ID[dof_Y];
				lm[2] = ss.Node(j).m_ID[dof_Z];
				lm[3] = ss.Node(j).m_ID[dof_RU];
				lm[4] = ss.Node(j).m_ID[dof_RV];
				lm[5] = ss.Node(j).m_ID[dof_RW];

				for (int k = 0; k < n; ++k)
				{
					vector<int>& id = ms.Node(en[k]).m_ID;
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
	}
	else
	{
		vector<int> lm;
		for (int j = 0; j < ss.Nodes(); ++j)
		{
			FESurfaceElement* pe = ss.m_data[j].m_pme;
			if (pe != 0)
			{
				FESurfaceElement& me = *pe;
				int* en = &me.m_lnode[0];

				int n = me.Nodes();
				lm.assign(3*(n+2), -1);

				lm[0] = ss.Node(j).m_ID[dof_X];
				lm[1] = ss.Node(j).m_ID[dof_Y];
				lm[2] = ss.Node(j).m_ID[dof_Z];

				for (int k = 0; k < n; ++k)
				{
					vector<int>& id = ms.Node(en[k]).m_ID;
					lm[3 * (k + 1)    ] = id[dof_X];
					lm[3 * (k + 1) + 1] = id[dof_Y];
					lm[3 * (k + 1) + 2] = id[dof_Z];
				}

				lm[3 * (n + 1)  ] = m_LM[3 * j   ];
				lm[3 * (n + 1)+1] = m_LM[3 * j +1];
				lm[3 * (n + 1)+2] = m_LM[3 * j +2];

				K.build_add(lm);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Interface activation
void FETiedInterface::Activate()
{
	// Don't forget to call base member!
	FEContactInterface::Activate();

	// project primary surface onto secondary surface
	ProjectSurface(ss, ms, m_breloc);
}

//-----------------------------------------------------------------------------
//! return number of equations to be allocated for Lagrange multipliers
int FETiedInterface::InitEquations(int neq)
{
	// make sure we want to use Lagrange Multiplier method
	if (m_laugon != 2) return 0;

	// allocate three equations per primary node
	int NN = ss.Nodes();

	m_LM.resize(3 * NN);
	for (int i = 0; i < 3 * NN; ++i) m_LM[i] = neq++;

	return 3 * NN;
}

//-----------------------------------------------------------------------------
//! Update tied interface data. This function re-evaluates the gaps between
//! the primary node and their projections onto the secondary surface.
//!
void FETiedInterface::Update()
{
	// get the mesh
	FEMesh& mesh = *ss.GetMesh();

	// loop over all primary nodes
	for (int i=0; i<ss.Nodes(); ++i)
	{
		FESurfaceElement* pme = ss.m_data[i].m_pme;
		if (pme)
		{
			// get the current primary nodal position
			vec3d rt = ss.Node(i).m_rt;

			// get the natural coordinates of the primary projection
			// onto the secondary element
			double r = ss.m_data[i].m_rs[0];
			double s = ss.m_data[i].m_rs[1];

			// get the nodal coordinates
			int ne = pme->Nodes();
			vec3d y[FEElement::MAX_NODES];
			for (int l=0; l<ne; ++l) y[l] = ms.Node( pme->m_lnode[l] ).m_rt;

			// calculate the primary node projection
			vec3d q = pme->eval(y, r, s);

			// calculate the secondary normal
			vec3d nu = ss.SurfaceNormal(*pme, r, s);

			// calculate the gap function
			// (taking possible offset into account)
			ss.m_data[i].m_vgap = (rt - q) - nu*ss.m_data[i].m_off;
			ss.m_data[i].m_gap = ss.m_data[i].m_vgap.norm();

			// calculate force
			ss.m_data[i].m_Tc = ss.m_data[i].m_Lm + ss.m_data[i].m_vgap*m_eps;
		}
	}
}

//-----------------------------------------------------------------------------
//! project surface

void FETiedInterface::ProjectSurface(FETiedContactSurface& ss, FETiedContactSurface& ms, bool bmove)
{
	// closest point projection method
	FEClosestPointProjection cpp(ms);
	cpp.SetTolerance(m_stol);
	cpp.HandleSpecialCases(m_bspecial);
	cpp.Init();

	// loop over all primary nodes
	for (int i=0; i<ss.Nodes(); ++i)
	{
		// get the next node
		FENode& node = ss.Node(i);
		ss.m_data[i].m_pme = nullptr;

		// get the nodal position of this primary node
		vec3d x = node.m_rt;

		// find the secondary element
		vec3d q; vec2d rs;
		FESurfaceElement* pme = cpp.Project(x, q, rs);
		if (pme)
		{
			// make sure we are within the max distance
			double D = (x - q).norm();
			if ((m_Dmax == 0.0) || (D <= m_Dmax))
			{
				// store the secondary element
				ss.m_data[i].m_pme = pme;

				// store the natural coordinates of the projection on the secondary element
				ss.m_data[i].m_rs = rs;

				// calculate the secondary normal
				vec3d nu = ms.SurfaceNormal(*pme, rs[0], rs[1]);

				// calculate gap
				ss.m_data[i].m_vgap = (x - q) - nu*ss.m_data[i].m_off;

				// move the node if necessary
				if (bmove && (ss.m_data[i].m_vgap.norm()>0))
				{
					node.m_r0 = node.m_rt = q + nu*ss.m_data[i].m_off;
					ss.m_data[i].m_vgap = vec3d(0,0,0);
				}

				// calculate force
				ss.m_data[i].m_Tc = ss.m_data[i].m_Lm + ss.m_data[i].m_vgap*m_eps;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact forces for a tied interface.

void FETiedInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
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

	// loop over all primary facets
	const int NE = ss.Elements();
	for (int j = 0; j < NE; ++j)
	{
		// get the primary element
		FESurfaceElement& sel = ss.Element(j);

		// get the element's LM vector
		ss.UnpackLM(sel, sLM);

		int nseln = sel.Nodes();

		double* w = sel.GaussWeights();

		// loop over primary element nodes (which are the integration points as well)
		for (int n = 0; n < nseln; ++n)
		{
			int m = sel.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a secondary element associated with it
			// TODO: is this a good way to test for an active constraint
			// The rigid wall criteria seems to work much better.
			if (ss.m_data[m].m_pme != 0)
			{
				// calculate jacobian and weight
				double Jw = ss.jac0(sel, n)*w[n];

				// get nodal contact force
				vec3d tc = ss.m_data[m].m_Lm;

				// add penalty contribution for penalty and aug lag method
				if (m_laugon != 2) tc += ss.m_data[m].m_vgap*m_eps;

				// get the secondary element
				FESurfaceElement& mel = *ss.m_data[m].m_pme;
				ms.UnpackLM(mel, mLM);
				int nmeln = mel.Nodes();

				// isoparametric coordinates of the projected primary node
				// onto the secondary element
				double r = ss.m_data[m].m_rs[0];
				double s = ss.m_data[m].m_rs[1];

				// get the secondary shape function values at this primary node
				mel.shape_fnc(N, r, s);

				// allocate "element" force vector
				if (m_laugon != 2) fe.resize(3 * (nmeln + 1));
				else fe.resize(3 * (nmeln + 2));

				// calculate contribution to force vector from nodes
				fe[0] = -Jw * tc.x;
				fe[1] = -Jw * tc.y;
				fe[2] = -Jw * tc.z;
				for (int l = 0; l < nmeln; ++l)
				{
					fe[3 * (l + 1)    ] = Jw * tc.x*N[l];
					fe[3 * (l + 1) + 1] = Jw * tc.y*N[l];
					fe[3 * (l + 1) + 2] = Jw * tc.z*N[l];
				}

				// setup lm vector
				if (m_laugon != 2) lm.resize(3 * (nmeln + 1));
				else lm.resize(3 * (nmeln + 2));

				// fill the lm array
				lm[0] = sLM[n * 3];
				lm[1] = sLM[n * 3 + 1];
				lm[2] = sLM[n * 3 + 2];
				for (int l = 0; l < nmeln; ++l)
				{
					lm[3 * (l + 1)    ] = mLM[l * 3];
					lm[3 * (l + 1) + 1] = mLM[l * 3 + 1];
					lm[3 * (l + 1) + 2] = mLM[l * 3 + 2];
				}

				if (m_laugon != 2) en.resize(nmeln + 1);
				else en.resize(nmeln + 2);

				// fill the en array
				en[0] = sel.m_node[n];
				for (int l = 0; l < nmeln; ++l) en[l + 1] = mel.m_node[l];

				if (m_laugon == 2)
				{
					// get the gap function
					vec3d g = ss.m_data[m].m_vgap;

					// add contribution from Lagrange multipliers
					fe[3 * (nmeln + 1)  ] = -Jw * g.x;
					fe[3 * (nmeln + 1)+1] = -Jw * g.y;
					fe[3 * (nmeln + 1)+2] = -Jw * g.z;

					// add the Lagrange multiplier equations to lm
					lm[3 * (nmeln + 1)  ] = m_LM[3 * m  ];
					lm[3 * (nmeln + 1)+1] = m_LM[3 * m+1];
					lm[3 * (nmeln + 1)+2] = m_LM[3 * m+2];

					// fill the en array
					en[nmeln + 1] = -1;
				}

				// assemble into global force vector
				R.Assemble(en, lm, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness matrix contribution.
void FETiedInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	vector<int> sLM, mLM, lm, en;
	const int MN = FEElement::MAX_NODES;
	FEElementMatrix ke;

	// shape functions
	double H[MN];

	// loop over all primary elements
	const int NE = ss.Elements();
	for (int i = 0; i < NE; ++i)
	{
		// get the next element
		FESurfaceElement& se = ss.Element(i);
		int nseln = se.Nodes();

		// get the element's LM vector
		ss.UnpackLM(se, sLM);

		double* w = se.GaussWeights();

		// loop over all integration points (that is nodes)
		for (int n = 0; n < nseln; ++n)
		{
			int m = se.m_lnode[n];

			// get the secondary element
			FESurfaceElement* pme = ss.m_data[m].m_pme;
			if (pme)
			{
				// get the secondary element
				FESurfaceElement& me = *pme;
				int nmeln = me.Nodes();
				ms.UnpackLM(me, mLM);

				// calculate jacobian
				double Jw = ss.jac0(se, n)*w[n];

				// primary node natural coordinates in secondary element
				double r = ss.m_data[m].m_rs[0];
				double s = ss.m_data[m].m_rs[1];

				// get the secondary shape function values at this primary node
				me.shape_fnc(H, r, s);

				if (m_laugon != 2)
				{
					// number of degrees of freedom
					int ndof = 3 * (1 + nmeln);

					// fill stiffness matrix
					ke.resize(ndof, ndof); ke.zero();
					ke[0][0] = Jw*m_eps;
					ke[1][1] = Jw*m_eps;
					ke[2][2] = Jw*m_eps;
					for (int k = 0; k < nmeln; ++k)
					{
						ke[0][3 + 3 * k    ] = -Jw*m_eps*H[k];
						ke[1][3 + 3 * k + 1] = -Jw*m_eps*H[k];
						ke[2][3 + 3 * k + 2] = -Jw*m_eps*H[k];

						ke[3 + 3 * k    ][0] = -Jw*m_eps*H[k];
						ke[3 + 3 * k + 1][1] = -Jw*m_eps*H[k];
						ke[3 + 3 * k + 2][2] = -Jw*m_eps*H[k];
					}
					for (int k = 0; k < nmeln; ++k)
						for (int l = 0; l < nmeln; ++l)
						{
							ke[3 + 3 * k    ][3 + 3 * l    ] = Jw*m_eps*H[k] * H[l];
							ke[3 + 3 * k + 1][3 + 3 * l + 1] = Jw*m_eps*H[k] * H[l];
							ke[3 + 3 * k + 2][3 + 3 * l + 2] = Jw*m_eps*H[k] * H[l];
						}

				}
				else 
				{
					// number of degrees of freedom
					int ndof = 3 * (2 + nmeln);

					// fill stiffness matrix
					ke.resize(ndof, ndof); ke.zero();

					int L = 3 * (nmeln + 1);
					ke[0][L  ] = ke[L  ][0] = Jw;
					ke[1][L+1] = ke[L+1][1] = Jw;
					ke[2][L+2] = ke[L+2][2] = Jw;
					for (int k = 0; k < nmeln; ++k)
					{
						ke[3 + 3*k    ][L    ] = -Jw*H[k];
						ke[3 + 3*k + 1][L + 1] = -Jw*H[k];
						ke[3 + 3*k + 2][L + 2] = -Jw*H[k];

						ke[L  ][3 + 3*k    ] = -Jw*H[k];
						ke[L+1][3 + 3*k + 1] = -Jw*H[k];
						ke[L+2][3 + 3*k + 2] = -Jw*H[k];
					}
				}

				// create lm array
				if (m_laugon != 2) lm.resize(3 * (1 + nmeln));
				else lm.resize(3 * (2 + nmeln));

				lm[0] = sLM[n * 3];
				lm[1] = sLM[n * 3 + 1];
				lm[2] = sLM[n * 3 + 2];
				for (int k = 0; k < nmeln; ++k)
				{
					lm[3 * (k + 1)] = mLM[k * 3];
					lm[3 * (k + 1) + 1] = mLM[k * 3 + 1];
					lm[3 * (k + 1) + 2] = mLM[k * 3 + 2];
				}

				if (m_laugon == 2)
				{
					lm[3 * (nmeln + 1)] = m_LM[3 * m];
					lm[3 * (nmeln + 1) + 1] = m_LM[3 * m + 1];
					lm[3 * (nmeln + 1) + 2] = m_LM[3 * m + 2];
				}

				// create the en array
				if (m_laugon != 2) en.resize(nmeln + 1);
				else en.resize(nmeln + 2);

				en[0] = se.m_node[n];
				for (int k = 0; k < nmeln; ++k) en[k + 1] = me.m_node[k];

				if (m_laugon == 2) en[nmeln + 1] = -1;

				// assemble stiffness matrix
				ke.SetNodes(en);
				ke.SetIndices(lm);
				LS.Assemble(ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Do an augmentation.
bool FETiedInterface::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != 1) return true;

	int i;

	// calculate initial norms
	double normL0 = 0;
	for (i=0; i<ss.Nodes(); ++i)
	{
		vec3d lm = ss.m_data[i].m_Lm;
		normL0 += lm*lm;
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i=0; i<ss.Nodes(); ++i)
	{
		vec3d lm = ss.m_data[i].m_Lm + ss.m_data[i].m_vgap*m_eps;

		normL1 += lm*lm;
		if (ss.m_data[i].m_pme != 0)
		{
			double g = ss.m_data[i].m_vgap.norm();
			normgc += g*g;
			++N;
		}
	}	
	if (N == 0) N=1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// check convergence of constraints
	feLog(" tied interface # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	feLog("    normal force : %15le %15le\n", pctn, m_atol);
	feLog("    gap function : %15le       ***\n", normgc);
		
	// check convergence
	bool bconv = true;
	if (pctn >= m_atol) bconv = false;
	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;

	if (bconv == false) 
	{
		for (i=0; i<ss.Nodes(); ++i)
		{
			// update Lagrange multipliers
			ss.m_data[i].m_Lm = ss.m_data[i].m_Lm + ss.m_data[i].m_vgap*m_eps;
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Serialize the data to the archive.
void FETiedInterface::Serialize(DumpStream &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	ms.Serialize(ar);
	ss.Serialize(ar);

	// serialize pointers
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			int NN = ss.Nodes();
			ar << NN;
			for (int i=0; i<NN; ++i)
			{
				FESurfaceElement* pe = ss.m_data[i].m_pme;
				if (pe) ar << pe->m_lid; else ar << -1;
			}
		}
		else
		{
			int NN, lid;
			ar >> NN;
			for (int i=0; i<NN; ++i)
			{
				ar >> lid;
				if (lid < 0) ss.m_data[i].m_pme = nullptr; else ss.m_data[i].m_pme = &ms.Element(lid);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Update Lagrange multipliers
void FETiedInterface::Update(vector<double>& ui)
{
	if (m_laugon == 2)
	{
		for (int i = 0; i < ss.Nodes(); ++i)
		{
			ss.m_data[i].m_Lm.x += ui[m_LM[3 * i    ]];
			ss.m_data[i].m_Lm.y += ui[m_LM[3 * i + 1]];
			ss.m_data[i].m_Lm.z += ui[m_LM[3 * i + 2]];
		}
	}
}
