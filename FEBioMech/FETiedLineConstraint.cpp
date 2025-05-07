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
#include "FETiedLineConstraint.h"
#include <FECore/FEClosestPointProjection.h>
#include <FECore/FELinearSystem.h>
#include <FECore/log.h>

FETiedLine::FETiedLine(FEModel* fem) : FEEdge(fem)
{

}

FEMaterialPoint* FETiedLine::CreateMaterialPoint()
{
	return new FELineMaterialPoint;
}

void FETiedLine::Update()
{

}

bool FETiedLine::Create(FESegmentSet& eset)
{
	return FEEdge::Create(eset, FE_LINE2NI); // TODO: Does this assume linear edges?
}

bool FETiedLine::Init()
{
	if (FEEdge::Init() == false) return false;
	int NN = Nodes();
	m_data.resize(NN);
	return true;
}

FETiedLine::Projection FETiedLine::ClosestProjection(const vec3d& x)
{
	FETiedLine::Projection P;
	int NE = Elements();
	double L2min = 0;
	const double eps = 1e-5;
	for (int i = 0; i < NE; ++i)
	{
		FELineElement& el = Element(i);
		vec3d a = Node(el.m_lnode[0]).m_rt;
		vec3d b = Node(el.m_lnode[1]).m_rt;
		vec3d e = (b - a);

		double D2 = e.norm2();
		if (D2 != 0)
		{
			double l = ((x-a)*e)/D2;
			if ((l >= -eps) && (l <= 1.0 + eps))
			{
				vec3d q = a + e * l;
				double L2 = (x - q).norm2();
				if ((P.pe == nullptr) || (L2 < L2min))
				{
					P.pe = &el;
					P.r = 2.0 * l - 1;
					P.q = q;
					L2min = L2;
				}
			}
		}
	}
	assert(P.pe);
	return P;
}

BEGIN_FECORE_CLASS(FETiedLineConstraint, FENLConstraint)
	ADD_PARAMETER(m_laugon  , "laugon"          )->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
	ADD_PARAMETER(m_atol    , "tolerance"       );
	ADD_PARAMETER(m_eps     , "penalty"         );
	ADD_PARAMETER(m_naugmin , "minaug"          );
	ADD_PARAMETER(m_naugmax , "maxaug"          );
	ADD_PARAMETER(m_stol    , "search_tolerance");
	ADD_PARAMETER(m_boffset , "offset_shells"   );
	ADD_PARAMETER(m_Dmax    , "max_distance"    );
	ADD_PARAMETER(m_bspecial, "special"         );
	ADD_PARAMETER(m_breloc  , "node_reloc"      );

	ADD_PROPERTY(pl, "primary")->AddFlag(FEProperty::Reference);
	ADD_PROPERTY(sl, "secondary")->AddFlag(FEProperty::Reference);

END_FECORE_CLASS();

FETiedLineConstraint::FETiedLineConstraint(FEModel* pfem) : FENLConstraint(pfem), pl(pfem), sl(pfem)
{
	static int count = 1;
	SetID(count++);

	// initial parameter values
	m_laugon = FECore::PENALTY_METHOD;
	m_atol = 0.01;
	m_eps = 1.0;
	m_stol = 0.0001;
	m_naugmin = 0;
	m_naugmax = 10;
	m_boffset = false;
	m_Dmax = 0.0;
	m_bspecial = true;
	m_breloc = false;
}

bool FETiedLineConstraint::Init()
{
	if (pl.Init() == false) return false;
	if (sl.Init() == false) return false;

	return true;
}

void FETiedLineConstraint::BuildMatrixProfile(FEGlobalMatrix& K)
{
	FEMesh& mesh = GetMesh();

	// get the DOFS
	const int dof_X = GetDOFIndex("x");
	const int dof_Y = GetDOFIndex("y");
	const int dof_Z = GetDOFIndex("z");
	const int dof_RU = GetDOFIndex("Ru");
	const int dof_RV = GetDOFIndex("Rv");
	const int dof_RW = GetDOFIndex("Rw");

	const int LMSIZE = 6 * (FEElement::MAX_NODES + 1);
	vector<int> lm(LMSIZE);

	for (int j = 0; j < pl.Nodes(); ++j)
	{
		FENode& nj = pl.Node(j);
		FELineElement* pe = pl.m_data[j].me;
		if (pe != 0)
		{
			FELineElement& me = *pe;
			int* en = &me.m_lnode[0];

			int n = me.Nodes();
			lm.assign(LMSIZE, -1);

			lm[0] = nj.m_ID[dof_X];
			lm[1] = nj.m_ID[dof_Y];
			lm[2] = nj.m_ID[dof_Z];
			lm[3] = nj.m_ID[dof_RU];
			lm[4] = nj.m_ID[dof_RV];
			lm[5] = nj.m_ID[dof_RW];

			for (int k = 0; k < n; ++k)
			{
				vector<int>& id = sl.Node(en[k]).m_ID;
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

void FETiedLineConstraint::Activate()
{
	// Don't forget to call base member!
	FENLConstraint::Activate();

	// project primary surface onto secondary surface
	ProjectLines(pl, sl);
}

void FETiedLineConstraint::Update()
{
	// get the mesh
	FEMesh& mesh = GetMesh();

	// loop over all primary nodes
	for (int i=0; i<pl.Nodes(); ++i)
	{
		FELineElement* pme = pl.m_data[i].me;
		if (pme)
		{
			// get the current primary nodal position
			vec3d rt = pl.Node(i).m_rt;

			// get the natural coordinates of the primary projection
			// onto the secondary element
			double r = pl.m_data[i].r;

			// get the nodal coordinates
			int ne = pme->Nodes();
			vec3d y[FEElement::MAX_NODES];
			for (int l=0; l<ne; ++l) y[l] = sl.Node( pme->m_lnode[l] ).m_rt;

			// calculate the primary node projection
			vec3d q;
			double H[FEElement::MAX_NODES];
			pme->shape(H, r);
			for (int n = 0; n < ne; ++n) q += y[n] * H[n];

			// calculate the gap function
			pl.m_data[i].vgap = (rt - q);

			// calculate force
			pl.m_data[i].Tc = pl.m_data[i].Lm + pl.m_data[i].vgap*m_eps;
		}
	}
}

void FETiedLineConstraint::ProjectLines(FETiedLine& pl, FETiedLine& sl)
{
	// let's count contact pairs
	int contacts = 0;

	// loop over all primary nodes
	for (int i=0; i< pl.Nodes(); ++i)
	{
		// get the next node
		FENode& node = pl.Node(i);
		pl.m_data[i].me = nullptr;

		// get the nodal position of this primary node
		vec3d x = node.m_rt;

		// find the secondary element
		FETiedLine::Projection P = sl.ClosestProjection(x);
		if (P.pe)
		{
			// make sure we are within the max distance
			double D = (x - P.q).norm();
			if ((m_Dmax == 0.0) || (D <= m_Dmax))
			{
				// store the secondary element
				pl.m_data[i].me = P.pe;

				// store the natural coordinates of the projection on the secondary element
				pl.m_data[i].r = P.r;

				// calculate vector gap
				pl.m_data[i].vgap = (x - P.q);

				// calculate force
				pl.m_data[i].Tc = pl.m_data[i].Lm + pl.m_data[i].vgap*m_eps;

				contacts++;
			}
		}
	}

	// if we found no contact pairs, let's report this since this is probably not the user's intention
	if (contacts == 0)
	{
		std::string name = GetName();
		feLogWarning("No contact pairs found for 1D slide line\"%s\".\nThis contact interface may not have any effect.", name.c_str());
	}
}

void FETiedLineConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
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
	const int NE = pl.Elements();
	for (int j = 0; j < NE; ++j)
	{
		// get the primary element
		FELineElement& sel = pl.Element(j);
		int nseln = sel.Nodes();

		// get the element's LM vector
		// TODO: This assumes dofs are indexed at (0,1,2)!
		sLM.resize(3 * nseln);
		for (int a = 0; a < nseln; ++a)
		{
			FENode& node = pl.Node(sel.m_lnode[a]);
			sLM[3 * a    ] = node.m_ID[0];
			sLM[3 * a + 1] = node.m_ID[1];
			sLM[3 * a + 2] = node.m_ID[2];
		}

		double* w = sel.GaussWeights();

		// loop over primary element nodes (which are the integration points as well)
		for (int n = 0; n < nseln; ++n)
		{
			int m = sel.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a secondary element associated with it
			// TODO: is this a good way to test for an active constraint
			// The rigid wall criteria seems to work much better.
			if (pl.m_data[m].me != 0)
			{
				// calculate jacobian and weight
				double Jw = 2*w[n]; // TODO: This assumes local jacobian is 2!

				// get nodal contact force
				vec3d tc = pl.m_data[m].Lm;

				// add penalty contribution for penalty and aug lag method
				tc += pl.m_data[m].vgap*m_eps;

				// get the secondary element
				FELineElement& mel = *pl.m_data[m].me;
				int nmeln = mel.Nodes();

				mLM.resize(3 * nmeln);
				for (int b = 0; b < nmeln; ++b)
				{
					FENode& node = sl.Node(mel.m_lnode[b]);
					mLM[3 * b    ] = node.m_ID[0];
					mLM[3 * b + 1] = node.m_ID[1];
					mLM[3 * b + 2] = node.m_ID[2];
				}

				// isoparametric coordinates of the projected primary node
				// onto the secondary element
				double r = pl.m_data[m].r;

				// get the secondary shape function values at this primary node
				mel.shape(N, r);

				// allocate "element" force vector
				fe.resize(3 * (nmeln + 1));

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
				lm.resize(3 * (nmeln + 1));

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

				en.resize(nmeln + 1);

				// fill the en array
				en[0] = sel.m_node[n];
				for (int l = 0; l < nmeln; ++l) en[l + 1] = mel.m_node[l];

				// assemble into global force vector
				R.Assemble(en, lm, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness matrix contribution.
void FETiedLineConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	vector<int> sLM, mLM, lm, en;
	const int MN = FEElement::MAX_NODES;
	FEElementMatrix ke;

	// shape functions
	double H[MN];

	// loop over all primary elements
	const int NE = pl.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FELineElement& se = pl.Element(i);
		int nseln = se.Nodes();

		// get the element's LM vector
		// TODO: This assumes dofs are indexed at (0,1,2)!
		sLM.resize(3 * nseln);
		for (int a = 0; a < nseln; ++a)
		{
			FENode& node = pl.Node(se.m_lnode[a]);
			sLM[3 * a    ] = node.m_ID[0];
			sLM[3 * a + 1] = node.m_ID[1];
			sLM[3 * a + 2] = node.m_ID[2];
		}

		double* w = se.GaussWeights();

		// loop over all integration points (that is nodes)
		for (int n = 0; n < nseln; ++n)
		{
			int m = se.m_lnode[n];

			// get the secondary element
			FELineElement* pme = pl.m_data[m].me;
			if (pme)
			{
				// get the secondary element
				FELineElement& me = *pme;
				int nmeln = me.Nodes();

				mLM.resize(3 * nmeln);
				for (int b = 0; b < nmeln; ++b)
				{
					FENode& node = sl.Node(me.m_lnode[b]);
					mLM[3 * b    ] = node.m_ID[0];
					mLM[3 * b + 1] = node.m_ID[1];
					mLM[3 * b + 2] = node.m_ID[2];
				}

				// calculate jacobian
				double Jw = 2*w[n]; // TODO: This assumes local jacobian is 2!

				// primary node natural coordinates in secondary element
				double r = pl.m_data[m].r;

				// get the secondary shape function values at this primary node
				me.shape(H, r);

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

				// create lm array
				lm.resize(3 * (1 + nmeln));

				lm[0] = sLM[n * 3];
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

bool FETiedLineConstraint::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != FECore::AUGLAG_METHOD) return true;

	// calculate initial norms
	double normL0 = 0;
	for (int i=0; i< pl.Nodes(); ++i)
	{
		vec3d lm = pl.m_data[i].Lm;
		normL0 += lm*lm;
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (int i=0; i<pl.Nodes(); ++i)
	{
		vec3d lm = pl.m_data[i].Lm + pl.m_data[i].vgap*m_eps;

		normL1 += lm*lm;
		if (pl.m_data[i].me != 0)
		{
			double g2 = pl.m_data[i].vgap.norm2();
			normgc += g2;
			++N;
		}
	}	
	if (N == 0) N=1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// check convergence of constraints
	const std::string& name = GetName();
	if (name.empty())
		feLog(" tied-line # %d\n", GetID());
	else
		feLog(" %s [tied-line]\n", name.c_str());

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
		for (int i=0; i<pl.Nodes(); ++i)
		{
			// update Lagrange multipliers
			pl.m_data[i].Lm = pl.m_data[i].Lm + pl.m_data[i].vgap*m_eps;
		}	
	}

	return bconv;
}

void FETiedLineConstraint::Serialize(DumpStream &ar)
{
	// store contact data
	FENLConstraint::Serialize(ar);

	// store contact surface data
	pl.Serialize(ar);
	sl.Serialize(ar);

	// serialize pointers
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			int NN = pl.Nodes();
			ar << NN;
			for (int i=0; i<NN; ++i)
			{
				FELineElement* pe = pl.m_data[i].me;
				if (pe) ar << pe->GetLocalID(); else ar << -1;
			}
		}
		else
		{
			int NN, lid;
			ar >> NN;
			for (int i=0; i<NN; ++i)
			{
				ar >> lid;
				if (lid < 0) pl.m_data[i].me = nullptr; 
				else pl.m_data[i].me = &sl.Element(lid);
			}
		}
	}
}
