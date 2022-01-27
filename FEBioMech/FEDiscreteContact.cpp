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
#include "FEDiscreteContact.h"
#include <FECore/FEModel.h>
#include <FECore/FEGlobalMatrix.h>
#include "FEContactInterface.h"
#include <FECore/FEClosestPointProjection.h>
#include <FECore/FELinearSystem.h>
#include <FECore/log.h>

FEDiscreteContactSurface::FEDiscreteContactSurface(FEModel* fem) : FEContactSurface(fem)
{

}

bool FEDiscreteContactSurface::Init()
{
	return FEContactSurface::Init();
}

BEGIN_FECORE_CLASS(FEDiscreteContact, FESurfaceConstraint)
	ADD_PARAMETER(m_blaugon, "laugon");
	ADD_PARAMETER(m_altol  , FE_RANGE_GREATER_OR_EQUAL(0.0), "altol");
	ADD_PARAMETER(m_gaptol , FE_RANGE_GREATER_OR_EQUAL(0.0), "gaptol");
	ADD_PARAMETER(m_penalty, FE_RANGE_GREATER_OR_EQUAL(0.0), "penalty");
	ADD_PARAMETER(m_naugmin, "minaug");
	ADD_PARAMETER(m_naugmax, "maxaug");
	ADD_PARAMETER(m_nsegup , "segup");
END_FECORE_CLASS();

FEDiscreteContact::FEDiscreteContact(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem)
{
	m_blaugon = false;
	m_altol = 0.01;
	m_penalty = 1.0;
	m_gaptol = 0.0;
	m_naugmin = 0;
	m_naugmax = 100;
	m_bfirst = true;
	m_nsegup = 0;
}

bool FEDiscreteContact::Init()
{
	return m_surf.Init();
}

void FEDiscreteContact::Activate()
{
	FENLConstraint::Activate();
	ProjectSurface(true);
}

void FEDiscreteContact::Update(const FETimeInfo& tp)
{
	bool bupdate = (m_bfirst || (m_nsegup == 0)? true : (tp.currentIteration <= m_nsegup));
	ProjectSurface(true);
	m_bfirst = false;
}

void FEDiscreteContact::SetDiscreteSet(FEDiscreteSet* pset)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	vector<int> tag(mesh.Nodes(), 0);
	m_Node.clear();
	int nsize = pset->size();
	for (int i=0; i<nsize; ++i)
	{
		const FEDiscreteSet::NodePair& delem = pset->Element(i);
		tag[delem.n0] = 1;
		tag[delem.n1] = 1;
	}
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		if (tag[i])
		{
			NODE node;
			node.nid = i;
			node.pe = 0;
			node.Lm = 0;
			node.gap = 0;
			node.nu = vec3d(0,0,0);
			node.q = vec3d(0,0,0);
			node.proj[0] = node.proj[1] = 0.0;

			m_Node.push_back(node);
		}
	}
}

void FEDiscreteContact::ProjectSurface(bool bsegup)
{
	FEClosestPointProjection cpp(m_surf);
	cpp.SetTolerance(0.01);
	cpp.SetSearchRadius(0.0);
	cpp.HandleSpecialCases(true);
	cpp.Init();

	// loop over all primary nodes
	FEMesh& mesh = *m_surf.GetMesh();
	for (int i=0; i<(int)m_Node.size(); ++i)
	{
		NODE& nodeData = m_Node[i];

		// get the node
		FENode& node = mesh.Node(nodeData.nid);

		// get the nodal position
		vec3d x = node.m_rt;

		// If the node is in contact, let's see if the node still is 
		// on the same secondary element
		if (nodeData.pe != 0)
		{
			FESurfaceElement& mel = *nodeData.pe;

			double r = nodeData.proj[0];
			double s = nodeData.proj[1];

			vec3d q = m_surf.ProjectToSurface(mel, x, r, s);
			nodeData.proj[0] = r;
			nodeData.proj[1] = s;
			nodeData.q = q;

			if (bsegup && (!m_surf.IsInsideElement(mel, r, s, 0.01)))
			{
				// see if the node might have moved to another secondary element
				vec2d rs(0,0);
				nodeData.pe = cpp.Project(x, q, rs);
				nodeData.proj[0] = rs.x();
				nodeData.proj[1] = rs.y();
				nodeData.q = q;
			}
		}
		else if (bsegup)
		{
			vec2d rs(0,0); vec3d q;
			nodeData.pe = cpp.Project(x, q, rs);
			nodeData.proj[0] = rs.x();
			nodeData.proj[1] = rs.y();
			nodeData.q = q;
		}

		// if we found a secondary element, update the gap and normal data
		if (nodeData.pe != 0)
		{
			FESurfaceElement& mel =  *nodeData.pe;

			double r = nodeData.proj[0];
			double s = nodeData.proj[1];

			// the primary normal is set to the secondary element normal
			nodeData.nu = m_surf.SurfaceNormal(mel, r, s);

			// calculate gap
			nodeData.gap = -(nodeData.nu*(x - nodeData.q));
		}
		else
		{
			// TODO: Is this a good criteria for out-of-contact?
			//		 perhaps this is not even necessary.
			// since the node is not in contact, we set the gap function 
			// and Lagrangian multiplier to zero
			nodeData.gap = 0;
			nodeData.Lm = 0;
		}
	}
}

void FEDiscreteContact::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	// the elements LM vectors
	vector<int> mLM;

	// loop over all primary nodes
	FEMesh& mesh = *m_surf.GetMesh();
	int nodes = (int) m_Node.size();
	for (int i=0; i<nodes; ++i)
	{
		NODE& nodeData = m_Node[i];

		FENode& node = mesh.Node(nodeData.nid);
		vector<int>& sLM = node.m_ID;

		FESurfaceElement* pe = nodeData.pe;

		// see if this node's constraint is active
		// that is, if it has a secondary element associated with it
		// TODO: is this a good way to test for an active constraint
		// The rigid wall criteria seems to work much better.
		if (pe != 0)
		{
			// This node is active and could lead to a non-zero
			// contact force.
			// get the secondary element
			FESurfaceElement& mel = *pe;
			m_surf.UnpackLM(mel, mLM);

			// calculate the degrees of freedom
			int nmeln = mel.Nodes();
			int ndof = 3*(nmeln+1);
			fe.resize(ndof);

			// calculate the nodal force
			ContactNodalForce(nodeData, mel, fe);

			// fill the lm array
			lm.resize(3*(nmeln+1));
			lm[0] = sLM[0];
			lm[1] = sLM[1];
			lm[2] = sLM[2];

			for (int l=0; l<nmeln; ++l)
			{
				lm[3*(l+1)  ] = mLM[l*3  ];
				lm[3*(l+1)+1] = mLM[l*3+1];
				lm[3*(l+1)+2] = mLM[l*3+2];
			}

			// fill the en array
			en.resize(nmeln+1);
			en[0] = nodeData.nid;
			for (int l=0; l<nmeln; ++l) en[l+1] = mel.m_node[l];

			// assemble into global force vector
			R.Assemble(en, lm, fe);
		}
	}
}

void FEDiscreteContact::ContactNodalForce(FEDiscreteContact::NODE& nodeData, FESurfaceElement& mel, vector<double>& fe)
{
	// max nr of secondary element nodes
	const int MAXMN = FEElement::MAX_NODES;

	// secondary element nodes
	vec3d rtm[MAXMN];

	// secondary shape function values at projection point
	double H[MAXMN];

	// contact forces
	double N[3*(MAXMN+1)];

	// get the mesh
	FEMesh& mesh = *m_surf.GetMesh();

	// gap function
	double gap = nodeData.gap;

	// penalty
	double eps = m_penalty;

	// get primary node normal force
	double Ln = nodeData.Lm;
	double tn = Ln + eps*gap;
	tn = MBRACKET(tn);

	// get the primary node normal
	vec3d nu = nodeData.nu;

	int nmeln = mel.Nodes();
	int ndof = 3*(1 + nmeln);

	// get the secondary element node positions
	for (int k=0; k<nmeln; ++k) rtm[k] = mesh.Node(mel.m_node[k]).m_rt;

	// isoparametric coordinates of the projected primary node
	// onto the secondary element
	double r = nodeData.proj[0];
	double s = nodeData.proj[1];

	// get the secondary shape function values at this primary node
	mel.shape_fnc(H, r, s);

	// calculate contact vectors for normal traction
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

void FEDiscreteContact::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEElementMatrix ke;

	const int MAXMN = FEElement::MAX_NODES;
	vector<int> lm(3*(MAXMN + 1));
	vector<int> en(MAXMN+1);

	vector<int> sLM;
	vector<int> mLM;

	// loop over all integration points (that is nodes)
	FEMesh& mesh = *m_surf.GetMesh();
	int nodes = (int) m_Node.size();
	for (int i=0; i<nodes; ++i)
	{
		NODE& nodeData = m_Node[i];

		vector<int>& sLM = mesh.Node(nodeData.nid).m_ID;

		// see if this node's constraint is active
		// that is, if it has a secondary element associated with it
		if (nodeData.pe != 0)
		{
			// get the secondary element
			FESurfaceElement& me = *nodeData.pe;

			// get the secondary surface element's LM array
			m_surf.UnpackLM(me, mLM);

			int nmeln = me.Nodes();
			int ndof = 3*(nmeln+1);

			// calculate the stiffness matrix
			ke.resize(ndof, ndof);
			ContactNodalStiffness(nodeData, me, ke);

			// fill the lm array
			lm[0] = sLM[0];
			lm[1] = sLM[1];
			lm[2] = sLM[2];
			for (int k=0; k<nmeln; ++k)
			{
				lm[3*(k+1)  ] = mLM[k*3  ];
				lm[3*(k+1)+1] = mLM[k*3+1];
				lm[3*(k+1)+2] = mLM[k*3+2];
			}

			// create the en array
			en.resize(nmeln+1);
			en[0] = nodeData.nid;
			for (int k=0; k<nmeln; ++k) en[k+1] = me.m_node[k];
						
			// assemble stiffness matrix
			ke.SetNodes(en);
			ke.SetIndices(lm);
			LS.Assemble(ke);
		}
	}
}

void FEDiscreteContact::ContactNodalStiffness(FEDiscreteContact::NODE& nodeData, FESurfaceElement& mel, matrix& ke)
{
	const int MAXMN = FEElement::MAX_NODES;

	vector<int> lm(3*(MAXMN+1));
	vector<int> en(MAXMN + 1);

	double H[MAXMN], Hr[MAXMN], Hs[MAXMN];
	double N[3*(MAXMN+1)], T1[3*(MAXMN+1)], T2[3*(MAXMN+1)];
	double N1[3*(MAXMN+1)], N2[3*(MAXMN+1)], D1[3*(MAXMN+1)], D2[3*(MAXMN+1)];
	double Nb1[3*(MAXMN+1)], Nb2[3*(MAXMN+1)];

	// get the mesh
	FEMesh& mesh = *m_surf.GetMesh();

	// nr of element nodes and degrees of freedom 
	int nmeln = mel.Nodes();
	int ndof = 3*(1 + nmeln);

	// penalty factor
	double eps = m_penalty;

	// nodal coordinates
	vec3d rt[MAXMN];
	for (int j=0; j<nmeln; ++j) rt[j] = mesh.Node(mel.m_node[j]).m_rt;

	// primary node natural coordinates in secondary element
	double r = nodeData.proj[0];
	double s = nodeData.proj[1];

	// primary gap
	double gap = nodeData.gap;

	// lagrange multiplier
	double Lm = nodeData.Lm;

	// get primary node normal force
	double tn = Lm + eps*gap;
	tn = MBRACKET(tn);

	// get the primary node normal
	vec3d nu = nodeData.nu;

	// get the secondary shape function values and the derivatives at this primary node
	mel.shape_fnc(H, r, s);
	mel.shape_deriv(Hr, Hs, r, s);

	// get the tangent vectors
	vec3d tau[2];
	m_surf.CoBaseVectors(mel, r, s, tau);

	// set up the N vector
	N[0] = nu.x;
	N[1] = nu.y;
	N[2] = nu.z;

	for (int k=0; k<nmeln; ++k) 
	{
		N[(k+1)*3  ] = -H[k]*nu.x;
		N[(k+1)*3+1] = -H[k]*nu.y;
		N[(k+1)*3+2] = -H[k]*nu.z;
	}
	
	// set up the Ti vectors
	T1[0] = tau[0].x; T2[0] = tau[1].x;
	T1[1] = tau[0].y; T2[1] = tau[1].y;
	T1[2] = tau[0].z; T2[2] = tau[1].z;

	for (int k=0; k<nmeln; ++k) 
	{
		T1[(k+1)*3  ] = -H[k]*tau[0].x;
		T1[(k+1)*3+1] = -H[k]*tau[0].y;
		T1[(k+1)*3+2] = -H[k]*tau[0].z;

		T2[(k+1)*3  ] = -H[k]*tau[1].x;
		T2[(k+1)*3+1] = -H[k]*tau[1].y;
		T2[(k+1)*3+2] = -H[k]*tau[1].z;
	}

	// set up the Ni vectors
	N1[0] = N2[0] = 0;
	N1[1] = N2[1] = 0;
	N1[2] = N2[2] = 0;

	for (int k=0; k<nmeln; ++k) 
	{
		N1[(k+1)*3  ] = -Hr[k]*nu.x;
		N1[(k+1)*3+1] = -Hr[k]*nu.y;
		N1[(k+1)*3+2] = -Hr[k]*nu.z;

		N2[(k+1)*3  ] = -Hs[k]*nu.x;
		N2[(k+1)*3+1] = -Hs[k]*nu.y;
		N2[(k+1)*3+2] = -Hs[k]*nu.z;
	}

	// calculate metric tensor
	mat2d M;
	M[0][0] = tau[0]*tau[0]; M[0][1] = tau[0]*tau[1]; 
	M[1][0] = tau[1]*tau[0]; M[1][1] = tau[1]*tau[1]; 

	// calculate reciprocal metric tensor
	mat2d Mi = M.inverse();

	// calculate curvature tensor
	double K[2][2] = {0};
	double Grr[FEElement::MAX_NODES];
	double Grs[FEElement::MAX_NODES];
	double Gss[FEElement::MAX_NODES];
	mel.shape_deriv2(Grr, Grs, Gss, r, s);
	for (int k=0; k<nmeln; ++k)
	{
		K[0][0] += (nu*rt[k])*Grr[k];
		K[0][1] += (nu*rt[k])*Grs[k];
		K[1][0] += (nu*rt[k])*Grs[k];
		K[1][1] += (nu*rt[k])*Gss[k];
	}

	// setup A matrix A = M + gK
	double A[2][2];
	A[0][0] = M[0][0] + gap*K[0][0];
	A[0][1] = M[0][1] + gap*K[0][1];
	A[1][0] = M[1][0] + gap*K[1][0];
	A[1][1] = M[1][1] + gap*K[1][1];

	// calculate determinant of A
	double detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];

	// setup Di vectors
	for (int k=0; k<ndof; ++k)
	{
		D1[k] = (1/detA)*(A[1][1]*(T1[k]+gap*N1[k]) - A[0][1]*(T2[k] + gap*N2[k]));
		D2[k] = (1/detA)*(A[0][0]*(T2[k]+gap*N2[k]) - A[0][1]*(T1[k] + gap*N1[k]));
	}

	// setup Nbi vectors
	for (int k=0; k<ndof; ++k)
	{
		Nb1[k] = N1[k] - K[0][1]*D2[k];
		Nb2[k] = N2[k] - K[0][1]*D1[k];
	}

	// --- N O R M A L   S T I F F N E S S ---
	double sum;
	for (int k=0; k<ndof; ++k)
		for (int l=0; l<ndof; ++l)
			{
				sum = 0;

				sum = Mi[0][0]*Nb1[k]*Nb1[l]+Mi[0][1]*(Nb1[k]*Nb2[l]+Nb2[k]*Nb1[l])+Mi[1][1]*Nb2[k]*Nb2[l];
				sum *= gap;
				sum -= D1[k]*N1[l]+D2[k]*N2[l]+N1[k]*D1[l]+N2[k]*D2[l];
				sum += K[0][1]*(D1[k]*D2[l]+D2[k]*D1[l]);
				sum *= tn;

				sum += eps*HEAVYSIDE(Lm+eps*gap)*N[k]*N[l];
	
				ke[k][l] = sum;
			}		
}


bool FEDiscreteContact::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (!m_blaugon) return true;

	bool bconv = true;
	mat2d Mi;

	// penalty factor
	double eps = m_penalty;

	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	double normL0 = 0;
	for (int i=0; i<(int) m_Node.size(); ++i)	normL0 += m_Node[i].Lm * m_Node[i].Lm;
	normL0 = sqrt(normL0);

	// --- c a l c u l a t e   c u r r e n t   n o r m s ---
	double normL1 = 0;	// force norm
	double normg1 = 0;	// gap norm
	int N = 0;
	for (int i=0; i<(int)m_Node.size(); ++i)
	{
		// update Lagrange multipliers
		double Ln = m_Node[i].Lm + eps*m_Node[i].gap;
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;

		if (m_Node[i].gap > 0)
		{
			normg1 += m_Node[i].gap*m_Node[i].gap;
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
	if (normg1 != 0) gnorm = fabs(normg1 - m_normg0)/normg1; else gnorm = fabs(normg1 - m_normg0);

	feLog(" discrete contact # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	feLog("    normal force : %15le", lnorm);
	if (m_altol > 0) feLog("%15le\n", m_altol); else feLog("       ***\n");
	feLog("    gap function : %15le", gnorm);
	if (m_gaptol > 0) feLog("%15le\n", m_gaptol); else feLog("       ***\n");

	// check convergence
	bconv = true;
	if ((m_altol > 0) && (lnorm > m_altol)) bconv = false;
	if ((m_gaptol > 0) && (gnorm > m_gaptol)) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if (m_naugmax <= naug) bconv = true;
		
	if (bconv == false)
	{
		// we did not converge so update multipliers
		for (int i=0; i<(int) m_Node.size(); ++i)
		{
			NODE& node = m_Node[i];
			// update Lagrange multipliers
			double Ln = node.Lm + eps*node.gap;
			node.Lm = MBRACKET(Ln);
		}	
	}

	// store the last gap norm
	m_normg0 = normg1;

	return bconv;
}

void FEDiscreteContact::BuildMatrixProfile(FEGlobalMatrix& K)
{
	// TODO: this is currently for max 6 nodes (hence 7=6+1)
	vector<int> lm(6*7);

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the DOFS
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");
	const int dof_RU = fem.GetDOFIndex("Ru");
	const int dof_RV = fem.GetDOFIndex("Rv");
	const int dof_RW = fem.GetDOFIndex("Rw");

	const int nodes = (int) m_Node.size();
	for (int i=0; i<nodes; ++i)
	{
		NODE& nodeData = m_Node[i];

		// get the FE node
		FENode& node = mesh.Node(nodeData.nid);

		// get the secondary surface element
		FESurfaceElement* pe = nodeData.pe;

		if (pe != 0)
		{
			FESurfaceElement& me = *pe;
			int* en = &me.m_node[0];

			// Note that we need to grab the rigid degrees of freedom as well
			// this is in case one of the nodes belongs to a rigid body.
			int n = me.Nodes();
			if (n == 3)
			{
				lm[6*(3+1)  ] = -1;lm[6*(3+2)  ] = -1;lm[6*(3+3)  ] = -1;
				lm[6*(3+1)+1] = -1;lm[6*(3+2)+1] = -1;lm[6*(3+3)+1] = -1;
				lm[6*(3+1)+2] = -1;lm[6*(3+2)+2] = -1;lm[6*(3+3)+2] = -1;
				lm[6*(3+1)+3] = -1;lm[6*(3+2)+3] = -1;lm[6*(3+3)+3] = -1;
				lm[6*(3+1)+4] = -1;lm[6*(3+2)+4] = -1;lm[6*(3+3)+4] = -1;
				lm[6*(3+1)+5] = -1;lm[6*(3+2)+5] = -1;lm[6*(3+3)+5] = -1;
			}
			if (n == 4)
			{
				lm[6*(4+1)  ] = -1;lm[6*(4+2)  ] = -1;
				lm[6*(4+1)+1] = -1;lm[6*(4+2)+1] = -1;
				lm[6*(4+1)+2] = -1;lm[6*(4+2)+2] = -1;
				lm[6*(4+1)+3] = -1;lm[6*(4+2)+3] = -1;
				lm[6*(4+1)+4] = -1;lm[6*(4+2)+4] = -1;
				lm[6*(4+1)+5] = -1;lm[6*(4+2)+5] = -1;
			}

			lm[0] = node.m_ID[dof_X];
			lm[1] = node.m_ID[dof_Y];
			lm[2] = node.m_ID[dof_Z];
			lm[3] = node.m_ID[dof_RU];
			lm[4] = node.m_ID[dof_RV];
			lm[5] = node.m_ID[dof_RW];

			for (int k=0; k<n; ++k)
			{
				vector<int>& id = mesh.Node(en[k]).m_ID;
				lm[6*(k+1)  ] = id[dof_X];
				lm[6*(k+1)+1] = id[dof_Y];
				lm[6*(k+1)+2] = id[dof_Z];
				lm[6*(k+1)+3] = id[dof_RU];
				lm[6*(k+1)+4] = id[dof_RV];
				lm[6*(k+1)+5] = id[dof_RW];
			}

			K.build_add(lm);
		}
	}
}

//=============================================================================
BEGIN_FECORE_CLASS(FEDiscreteContact2, FESurfaceConstraint);
END_FECORE_CLASS();

FEDiscreteContact2::FEDiscreteContact2(FEModel* fem) : FESurfaceConstraint(fem), m_surf(fem)
{
	m_dom = 0;	
}

bool FEDiscreteContact2::Init()
{
	// let's make sure we have a discrete domain set
	if (m_dom == 0) return false;

	// initialize node data
	int NN = m_dom->Nodes();
	m_nodeData.resize(NN);
	for (int i=0; i<NN; ++i)
	{
		NODE& nd = m_nodeData[i];
		nd.node = i;
		nd.pe = 0;
	}

	// initialize surface
	if (m_surf.Init() == false) return false;

	// let's do base class initiialization
	return FENLConstraint::Init();
}

void FEDiscreteContact2::Activate()
{
	FENLConstraint::Activate();
	ProjectNodes();
}

void FEDiscreteContact2::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	int NN = m_dom->Nodes();

	FEMesh& mesh = *m_surf.GetMesh();

	// skip first and last
	for (int n=1; n<NN-1; ++n)
	{
		if (m_dom->IsAnchored(n))
		{
			NODE& nodeData = m_nodeData[n];
			assert(nodeData.pe);
			FESurfaceElement& el = *nodeData.pe;

			double r = nodeData.proj[0];
			double s = nodeData.proj[1];

			double H[FEElement::MAX_NODES];
			int neln = el.Nodes();

			el.shape_fnc(H, r, s);

			vec3d F = m_dom->NodalForce(n);
			vec3d nu = nodeData.nu;
			if (F*nu < 0.0)
			{
				vector<double> fe(3*neln, 0.0);
				for (int i=0; i<neln; ++i)
				{
					fe[3*i  ] = H[i]*F.x;
					fe[3*i+1] = H[i]*F.y;
					fe[3*i+2] = H[i]*F.z;
				}

				vector<int> lm(3*neln, -1);
				for (int i=0; i<neln; ++i)
				{
					FENode& nodei = mesh.Node(el.m_node[i]);
					lm[3*i  ] = nodei.m_ID[0];
					lm[3*i+1] = nodei.m_ID[1];
					lm[3*i+2] = nodei.m_ID[2];
				}

				R.Assemble(el.m_node, lm, fe);
			}
		}
	}
}

void FEDiscreteContact2::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
}

void FEDiscreteContact2::BuildMatrixProfile(FEGlobalMatrix& M)
{
}

void FEDiscreteContact2::Update(const FETimeInfo& tp)
{
	ProjectNodes();
}

void FEDiscreteContact2::ProjectNodes()
{
	// setup closest point projection
	FEClosestPointProjection cpp(m_surf);
	cpp.SetTolerance(0.01);
	cpp.SetSearchRadius(0.0);
	cpp.HandleSpecialCases(true);
	cpp.Init();

	// number of nodes
	int NN = m_dom->Nodes();

	double L = m_dom->InitialLength();

	bool bdone = false;

	m_dom->UpdateNodes();

	const double snap_tol = 1e-6;
	const int maxIter = 1000;
	int iter = 0;
	while ((bdone == false) && (iter < maxIter))
	{
		bdone = true;
		iter++;

		// skip first and last
		for (int i = 1; i<NN - 1; ++i)
		{
			FENode& nd = m_dom->Node(i);
			NODE& nodeData = m_nodeData[i];
				
			// get nodal position
			vec3d x = nd.m_rt;

			// If the node is in contact, let's see if it remains in contact
			if (nodeData.pe != 0)
			{
				// TODO: project the node on the surface to make sure we are using the updated projection position
				FESurfaceElement& mel = *nodeData.pe;
				double r = nodeData.proj[0];
				double s = nodeData.proj[1];

				// find the position of neighbors
				vec3d ra = m_dom->Node(i - 1).m_rt;
				vec3d rb = m_dom->Node(i + 1).m_rt;

				// evaluate center
				vec3d c = (ra + rb)*0.5;

				// try the same element first
				FESurfaceElement* pe = nodeData.pe;
				vec3d q = m_surf.ProjectToSurface(mel, c, r, s);
				if (m_surf.IsInsideElement(mel, r, s, 0.05) == false)
				{
					// project onto surface
					vec2d rs(r,s);
					pe = cpp.Project(c, q, rs);
					r = rs.x();
					s = rs.y();
				}

				if (pe)
				{
					vec3d nu = m_surf.SurfaceNormal(*pe, r, s);
					double gap = nu*(c - q);
					if (gap < 0.0)
					{
						nodeData.pe = pe;
						nodeData.proj[0] = r;
						nodeData.proj[1] = s;
						nodeData.nu = nu;
						nodeData.q = q;

						// see if the node moved significantly
						double d = (q - x).norm();
						if (d / L > 1e-2) bdone = false;
					}
					else if (gap > snap_tol)
					{
						// projection failed, release node and redo
						m_dom->AnchorNode(i, false);
						nodeData.pe = 0;
						bdone = false;
					}
					else if ((x - c).norm2() > snap_tol)
					{
						// projection failed, release node and redo
						m_dom->AnchorNode(i, false);
						nodeData.pe = 0;
						bdone = false;
					}
				}
				else
				{
					// projection faild, release node and redo
					m_dom->AnchorNode(i, false);
					nodeData.pe = 0;
					bdone = false;
				}
			}
			else
			{
				// see it the node establishes contact
				vec2d rs(0, 0); vec3d q;
				FESurfaceElement* pe = cpp.Project(x, q, rs);

				// if it does, update nodal data
				// and position the node
				if (pe)
				{
					vec3d nu = m_surf.SurfaceNormal(*pe, rs.x(), rs.y());
					if (nu*(x - q) < 0)
					{
						nodeData.pe = pe;
						nodeData.proj[0] = rs.x();
						nodeData.proj[1] = rs.y();
						nodeData.nu = nu;
						nodeData.q = q;

						m_dom->AnchorNode(i, true);
						bdone = false;
					}
				}
			}
		}

		// update anchor positions
		for (int i=1; i<NN-1; ++i)
		{
			NODE& nodeData = m_nodeData[i];
			if (nodeData.pe != 0) m_dom->SetNodePosition(i, nodeData.q);
		}

		// tell the domain to update the positions of the free nodes
		m_dom->UpdateNodes();
	}

//	assert(iter < maxIter);
	feLog("iterations = %d\n", iter);
}
