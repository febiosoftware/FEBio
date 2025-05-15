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
#include "FESlideLineConstraint.h"
#include <FECore/FEShellDomain.h>
#include "FECore/FEClosestPointProjection.h"
#include "FECore/FEModel.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FEAnalysis.h>

FESlideLine::FESlidingPoint::FESlidingPoint()
{
	pme = nullptr;
	gap = 0;
	nu = vec3d(0, 0, 0);
	r = 0;
	Lm = 0.0;
	Ln = 0.0;
}

void FESlideLine::FESlidingPoint::Serialize(DumpStream& ar)
{
	ar & gap & nu;
	ar & r;
	ar & Lm & Ln;
}

void FESlideLine::FESlidingPoint::Init()
{
	pme = nullptr;
	gap = 0;
	nu = vec3d(0, 0, 0);
	r = 0;
	Lm = Ln = 0.0;
}

FESlideLine::FESlideLine(FEModel* pfem) : FEEdge(pfem) {}

bool FESlideLine::Init()
{
	// always intialize base class first!
	if (FEEdge::Init() == false) return false;

	// get the number of nodes
	int nn = Nodes();

	// allocate integration point data
	m_data.resize(nn);
	for (int i = 0; i < nn; ++i)
	{
		FESlidingPoint& d = m_data[i];
		d.Init();
	}

	return true;
}

FEMaterialPoint* FESlideLine::CreateMaterialPoint()
{
	return new FELineMaterialPoint;
}

FESlideLine::Projection FESlideLine::ClosestProjection(const vec3d& x)
{
	FESlideLine::Projection P;
	int NE = Elements();
	double L2min = 0;
	for (int i = 0; i < NE; ++i)
	{
		FELineElement& el = Element(i);
		vec3d a = Node(el.m_lnode[0]).m_rt;
		vec3d b = Node(el.m_lnode[1]).m_rt;
		vec3d e = (b - a);

		double D2 = e.norm2();
		if (D2 != 0)
		{
			double l = ((x - a) * e) / D2;
			if (l < 0) l = 0;
			else if (l > 1) l = 1;

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
	assert(P.pe);
	return P;
}

void FESlideLine::Serialize(DumpStream& ar)
{
	// serialize base class
	FEEdge::Serialize(ar);

	// serialize data
//	ar & m_data;
}

BEGIN_FECORE_CLASS(FESlideLineConstraint, FENLConstraint)
	ADD_PARAMETER(m_laugon       , "laugon"       )->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
	ADD_PARAMETER(m_atol         , "tolerance"    );
	ADD_PARAMETER(m_eps          , "penalty"      );
	ADD_PARAMETER(m_gtol         , "gaptol"       );
	ADD_PARAMETER(m_naugmin      , "minaug"       );
	ADD_PARAMETER(m_naugmax      , "maxaug"       );
	ADD_PARAMETER(m_nsegup       , "seg_up"       );

	ADD_PROPERTY(pl, "primary")->AddFlag(FEProperty::Reference);
	ADD_PROPERTY(sl, "secondary")->AddFlag(FEProperty::Reference);

END_FECORE_CLASS();

//! build the matrix profile for use in the stiffness matrix
void FESlideLineConstraint::BuildMatrixProfile(FEGlobalMatrix& K)
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
		FELineElement* pe = pl.m_data[j].pme;
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

FESlideLineConstraint::FESlideLineConstraint(FEModel* pfem) : FENLConstraint(pfem), pl(pfem), sl(pfem)
{
	static int count = 1;
	SetID(count++);

	m_laugon = 0; // penalty method by default
	m_naugmin = 0;
	m_naugmax = 10;

	m_gtol = 0;

	m_nsegup = 0;	// always do segment updates
};

bool FESlideLineConstraint::Init()
{
	// set data
	m_bfirst = true;
	m_normg0 = 0.0;

	// create the surfaces
	if (pl.Init() == false) return false;
	if (sl.Init() == false) return false;

	return true;
}

void FESlideLineConstraint::Activate()
{
	// don't forget to call the base class
	FENLConstraint::Activate();

	// project primary surface onto secondary surface
	ProjectPoints(true);
}

void FESlideLineConstraint::ProjectPoints(bool bupseg)
{
	// loop over all primary nodes
	for (int i=0; i<pl.Nodes(); ++i)
	{
		FENode& node = pl.Node(i);
		vec3d x = node.m_rt;

		FESlideLine::FESlidingPoint& pti = pl.m_data[i];

		FESlideLine::Projection P = sl.ClosestProjection(x);
		if (P.pe)
		{
			FELineElement* pme = pti.pme;

			if ((P.pe == pme) || bupseg)
			{
				pti.pme = P.pe;
				pti.r = P.r;
				// calculate gap
				vec3d e = P.q - x;
				pti.gap = e.unit();
				pti.nu = e;
			}
		}
		else
		{
			// point has left contact
			pti.pme = nullptr;
			pti.gap = 0;
			pti.Ln = pti.Lm = 0;
		}
	}
}

//-----------------------------------------------------------------------------
//! updates sliding interface data
//! niter is the number of Newton iterations.
void FESlideLineConstraint::Update()
{
	int niter = GetFEModel()->GetCurrentStep()->GetFESolver()->m_niter;

	// should we do a segment update or not?
	// TODO: check what happens when m_nsegup == -1 and m_npass = 2;
	// We have to make sure that in this case, both surfaces get at least
	// one pass!
	bool bupdate = (m_bfirst || (m_nsegup == 0)? true : (niter <= m_nsegup));

	ProjectPoints(bupdate);

	// Update the net contact pressures
	UpdateContactPressures();

	// set the first-entry-flag to false
	m_bfirst = false;
}

void FESlideLineConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
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
	double w[MN];
	double detJ[MN];

	// loop over all primary surface facets
	int ne = pl.Elements();
	for (int j=0; j<ne; ++j)
	{
		// get the next element
		FELineElement& sel = pl.Element(j);
		int nseln = sel.Nodes();

		vec3d a = pl.Node(sel.m_lnode[0]).m_r0;
		vec3d b = pl.Node(sel.m_lnode[1]).m_r0;
		double L = (b - a).norm();

		// get the element's LM array
		// TODO: This assumes dofs are indexed at (0,1,2)!
		sLM.resize(3 * nseln);
		for (int a = 0; a < nseln; ++a)
		{
			FENode& node = pl.Node(sel.m_lnode[a]);
			sLM[3 * a] = node.m_ID[0];
			sLM[3 * a + 1] = node.m_ID[1];
			sLM[3 * a + 2] = node.m_ID[2];
		}

		// we calculate all the metrics we need before we
		// calculate the nodal forces
		for (int n=0; n<nseln; ++n)
		{
			// jacobians
			detJ[n] = L/2; // TODO: This assumes local jacobian is 2!

			// integration weights
			w[n] = 1;// sel.GaussWeights()[n];
		}

		// loop over primary surface element nodes (which are the integration points as well)
		// and calculate the contact nodal force
		for (int n=0; n<nseln; ++n)
		{
			// get the local node number
			int m = sel.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has an element associated with it
			// TODO: is this a good way to test for an active constraint
			// The rigid wall criteria seems to work much better.
			if (pl.m_data[m].pme != 0)
			{
				// This node is active and could lead to a non-zero
				// contact force.
				// get the secondary surface element
				FELineElement& mel = *pl.m_data[m].pme;
				int nmeln = mel.Nodes();

				mLM.resize(3 * nmeln);
				for (int b = 0; b < nmeln; ++b)
				{
					FENode& node = sl.Node(mel.m_lnode[b]);
					mLM[3 * b] = node.m_ID[0];
					mLM[3 * b + 1] = node.m_ID[1];
					mLM[3 * b + 2] = node.m_ID[2];
				}

				// calculate the nodal force
				int ndof = 3 * (nmeln + 1);
				fe.resize(ndof);
				ContactNodalForce(m, mel, fe);

				// multiply force with weights
				for (int l=0; l<ndof; ++l) fe[l] *= detJ[n]*w[n];
					
				// fill the lm array
				lm.resize(3*(nmeln+1));
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				for (int l=0; l<nmeln; ++l)
				{
					lm[3*(l+1)  ] = mLM[l*3  ];
					lm[3*(l+1)+1] = mLM[l*3+1];
					lm[3*(l+1)+2] = mLM[l*3+2];
				}

				// fill the en array
				en.resize(nmeln+1);
				en[0] = sel.m_node[n];
				for (int l=0; l<nmeln; ++l) en[l+1] = mel.m_node[l];

				// assemble into global force vector
				R.Assemble(en, lm, fe);
			}
		}
	}
}

//! Calculates the contact force on a node.
//! \param[in] m local node number
//! \param[out] fe force vector
void FESlideLineConstraint::ContactNodalForce(int m, FELineElement& mel, vector<double>& fe)
{
	// max nr of element nodes
	const int MAXMN = FEElement::MAX_NODES;

	// get primary node force
	double gap = pl.m_data[m].gap;
	double eps = Penalty();
	double Ln = pl.m_data[m].Lm;
	double tn = Ln + eps*gap;

	// get the primary node normal
	vec3d& nu = pl.m_data[m].nu;

	int nmeln = mel.Nodes();
	int ndof = 3*(1 + nmeln);

	// isoparametric coordinates of the projected node
	// onto the secondary surface element
	double r = pl.m_data[m].r;

	// get the secondary surface element shape function values at this node
	double H[MAXMN];
	mel.shape(H, r);

	// --- N O R M A L   T R A C T I O N ---

	// calculate contact vectors for normal traction
	double N[3 * (MAXMN + 1)];
	N[0] = nu.x;
	N[1] = nu.y;
	N[2] = nu.z;
	for (int l=0; l<nmeln; ++l)
	{
		N[3*(l+1)  ] = -H[l]*nu.x;
		N[3*(l+1)+1] = -H[l]*nu.y;
		N[3*(l+1)+2] = -H[l]*nu.z;
	}

	// calculate force vector
	for (int l=0; l<ndof; ++l) fe[l] = tn*N[l];
}

void FESlideLineConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEElementMatrix ke;

	const int MAXMN = FEElement::MAX_NODES;
	vector<int> lm(3*(MAXMN + 1));
	vector<int> en(MAXMN+1);

	double w[MAXMN];
	double detJ[MAXMN];

	vector<int> sLM;
	vector<int> mLM;

	// loop over all primary elements
	int ne = pl.Elements();
	for (int j=0; j<ne; ++j)
	{
		// unpack the next element
		FELineElement& se = pl.Element(j);
		int nseln = se.Nodes();

		vec3d a = pl.Node(se.m_lnode[0]).m_r0;
		vec3d b = pl.Node(se.m_lnode[1]).m_r0;
		double L = (b - a).norm();

		// get the element's LM array
		// TODO: This assumes dofs are indexed at (0,1,2)!
		sLM.resize(3 * nseln);
		for (int a = 0; a < nseln; ++a)
		{
			FENode& node = pl.Node(se.m_lnode[a]);
			sLM[3 * a] = node.m_ID[0];
			sLM[3 * a + 1] = node.m_ID[1];
			sLM[3 * a + 2] = node.m_ID[2];
		}

		// get all the metrics we need 
		for (int n=0; n<nseln; ++n)
		{
			detJ[n] = L/2;
			w[n] = 1;
		}

		// loop over all integration points (that is nodes)
		for (int n=0; n<nseln; ++n)
		{
			int m = se.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has an element associated with it
			if (pl.m_data[m].pme != 0)
			{
				// get the secondary surface element
				FELineElement& me = *pl.m_data[m].pme;

				int nmeln = me.Nodes();
				int ndof = 3 * (nmeln + 1);

				// get the secondary surface element's LM array
				mLM.resize(3 * nmeln);
				for (int b = 0; b < nmeln; ++b)
				{
					FENode& node = sl.Node(me.m_lnode[b]);
					mLM[3 * b] = node.m_ID[0];
					mLM[3 * b + 1] = node.m_ID[1];
					mLM[3 * b + 2] = node.m_ID[2];
				}

				// calculate the stiffness matrix
				ke.resize(ndof, ndof);
				ContactNodalStiffness(m, me, ke);

				// muliply with weights
				for (int k=0; k<ndof; ++k)
					for (int l=0; l<ndof; ++l) ke[k][l] *= detJ[n]*w[n];

				// fill the lm array
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				for (int k=0; k<nmeln; ++k)
				{
					lm[3*(k+1)  ] = mLM[k*3  ];
					lm[3*(k+1)+1] = mLM[k*3+1];
					lm[3*(k+1)+2] = mLM[k*3+2];
				}

				// create the en array
				en.resize(nmeln+1);
				en[0] = se.m_node[n];
				for (int k=0; k<nmeln; ++k) en[k+1] = me.m_node[k];
						
				// assemble stiffness matrix
				ke.SetNodes(en);
				ke.SetIndices(lm);
				LS.Assemble(ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlideLineConstraint::ContactNodalStiffness(int m, FELineElement& mel, matrix& ke)
{
	const int MAXMN = FEElement::MAX_NODES;

	vector<int> lm(3*(MAXMN+1));
	vector<int> en(MAXMN + 1);

	double H[MAXMN], Hr[MAXMN];

	double N[3][3 * (MAXMN + 1)] = { 0 };
	double AN[3][3 * (MAXMN + 1)] = { 0 };

	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// nr of element nodes and degrees of freedom 
	int nmeln = mel.Nodes();
	int ndof = 3*(1 + nmeln);

	// penalty factor
	double eps = Penalty();

	// nodal coordinates
	vec3d rt[MAXMN];
	for (int j=0; j<nmeln; ++j) rt[j] = mesh.Node(mel.m_node[j]).m_rt;

	// node natural coordinates in secondary surface element
	double r = pl.m_data[m].r;

	// gap
	double gap = pl.m_data[m].gap;

	// lagrange multiplier
	double Lm = pl.m_data[m].Lm;

	// get node normal force
	double tn = Lm + eps*gap;

	// get the node normal
	vec3d& nu = pl.m_data[m].nu;

	// get the secondary surface element shape function values and the derivatives at this node
	mel.shape(H, r);
	mel.shape_deriv(Hr, r);

	// set up the N vector
	N[0][0] = 1;
	N[1][1] = 1;
	N[2][2] = 1;
	for (int k=0; k<nmeln; ++k)
	{
		N[0][(k+1)*3  ] = -H[k];
		N[1][(k+1)*3+1] = -H[k];
		N[2][(k+1)*3+2] = -H[k];
	}

	// get the tangent vector
	vec3d tau(0, 0, 0);
	for (int i = 0; i < nmeln; ++i)
	{
		tau += mesh.Node(mel.m_node[i]).m_rt * Hr[i];
	}

	// calculate curvature tensor
	double K = 0;
	double Grr[FEElement::MAX_NODES];
	mel.shape_deriv2(Grr, r);
	for (int k=0; k<nmeln; ++k)
	{
		K += (nu*rt[k])*Grr[k];
	}

	// metric
	double t = tau * tau;

	// setup A matrix A = M + gK
	double A = t - gap*K;

	mat3d TxT = (tau & tau);

	mat3d B = mat3dd(1.0) - TxT/A;

	mat3d M = B.transpose() * B;

	for (int i=0; i<3; ++i)
		for (int j = 0; j < ndof; ++j)
		{
			AN[i][j] = 0;
			for (int k=0; k<3; ++k)
				AN[i][j] += M[i][k] * N[k][j];
		}

	for (int k=0; k<ndof; ++k)
		for (int l=0; l<ndof; ++l)
			{
				ke[k][l] = 0;
				for (int i=0; i<3; ++i)
					ke[k][l] += N[i][k] * AN[i][l];

				ke[l][k] *= eps;
			}
}

bool FESlideLineConstraint::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon != FECore::AUGLAG_METHOD) return true;

	double Ln;
	bool bconv = true;
	mat2d Mi;

	// penalty factor
	double eps = Penalty();

	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0;
	for (int i=0; i<pl.Nodes(); ++i) normL0 += pl.m_data[i].Lm*pl.m_data[i].Lm;
	normL0 = sqrt(normL0);

	// --- c a l c u l a t e   c u r r e n t   n o r m s ---
	// a. normal component
	double normL1 = 0;	// force norm
	double normg1 = 0;	// gap norm
	int N = 0;
	for (int i=0; i<pl.Nodes(); ++i)
	{
		// update Lagrange multipliers
		Ln = pl.m_data[i].Lm + eps*pl.m_data[i].gap;
		normL1 += Ln*Ln;
		if (pl.m_data[i].gap > 0)
		{
			normg1 += pl.m_data[i].gap*pl.m_data[i].gap;
			++N;
		}
	}	
	if (N == 0) N=1;

	normL1 = sqrt(normL1);
	normg1 = sqrt(normg1 / N);

	if (naug == 0) m_normg0 = 0;

	// calculate and print convergence norms
	double lnorm = 0, gnorm = 0;
	if (normL1 != 0) lnorm = fabs(normL1 - normL0)/normL1; else lnorm = fabs(normL1 - normL0);
//	if (normg1 != 0) gnorm = fabs(normg1 - m_normg0)/normg1; else gnorm = fabs(normg1 - m_normg0);
	gnorm = fabs(normg1 - m_normg0);

	feLog(" slide line # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	feLog("    normal force : %15le", lnorm);
	if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
	feLog("    gap function : %15le", gnorm);
	if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");

	// check convergence
	bconv = true;
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_gtol > 0) && (gnorm > m_gtol)) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if (m_naugmax <= naug) bconv = true;
		
	if (bconv == false)
	{
		// we did not converge so update multipliers
		for (int i=0; i<pl.Nodes(); ++i)
		{
			// update Lagrange multipliers
			Ln = pl.m_data[i].Lm + eps*pl.m_data[i].gap;
			pl.m_data[i].Lm = Ln;
		}
	}

	// store the last gap norm
	m_normg0 = normg1;

	return bconv;
}

void FESlideLineConstraint::UpdateContactPressures()
{
	double eps = Penalty();
	// loop over all nodes of the primary line
	for (int n=0; n<pl.Nodes(); ++n)
	{
		// get the normal tractions at the integration points
		double gap = pl.m_data[n].gap;
		pl.m_data[n].Ln = pl.m_data[n].Lm + eps*gap;
	}
}

void FESlideLineConstraint::Serialize(DumpStream& ar)
{
	// store contact data
	FENLConstraint::Serialize(ar);

	// store contact surface data
	pl.Serialize(ar);
	sl.Serialize(ar);

	// restore element pointers
	SerializePointers(ar);

	ar & m_bfirst & m_normg0;
}

void FESlideLineConstraint::SerializePointers(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		for (int i = 0; i < pl.m_data.size(); i++)
		{
			FELineElement* pe = pl.m_data[i].pme;
			int eid = (pe ? pe->GetLocalID() : -1);
			ar << eid;
		}
	}
	else
	{
		for (int i = 0; i < pl.m_data.size(); i++)
		{
			int eid = -1;
			ar >> eid;
			if (eid >= 0)
			{
				FELineElement* pe = &sl.Element(eid);
				pl.m_data[i].pme = pe;
			}
			else pl.m_data[i].pme = nullptr;
		}
	}
}
