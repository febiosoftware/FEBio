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
#include "FEDeformableSpringDomain.h"
#include <FECore/FEModel.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FELinearSystem.h>
#include "FEBioMech.h"

BEGIN_FECORE_CLASS(FEDeformableSpringDomain, FEDiscreteDomain)
	ADD_PARAMETER(m_kbend, "k_bend");
	ADD_PARAMETER(m_kstab, "k_stab");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEDeformableSpringDomain::FEDeformableSpringDomain(FEModel* pfem) : FEDiscreteDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem), m_dofR(pfem), m_dof(pfem)
{
	m_pMat  =   0;
	m_kbend = 0.0;
	m_kstab = 0.0;

	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
	}
}

//-----------------------------------------------------------------------------
//! get the total dofs
const FEDofList& FEDeformableSpringDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEDeformableSpringDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FESpringMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
void FEDeformableSpringDomain::BuildMatrixProfile(FEGlobalMatrix& K)
{
	// we connect each node to its two neighbors
	int NN = Nodes();
	for (int i=1; i<NN-1; ++i)
	{
		vector<int> lm(3*6, -1);
		for (int j=0; j<3; ++j)
		{
			int n = i-1+j;
			vector<int>& id = Node(n).m_ID;

			// first the displacement dofs
			lm[6 * j    ] = id[m_dofU[0]];
			lm[6 * j + 1] = id[m_dofU[1]];
			lm[6 * j + 2] = id[m_dofU[2]];

			// rigid rotational dofs
			lm[6 * j + 3] = id[m_dofR[0]];
			lm[6 * j + 4] = id[m_dofR[1]];
			lm[6 * j + 5] = id[m_dofR[2]];
		}

		K.build_add(lm);
	}
}

//-----------------------------------------------------------------------------
void FEDeformableSpringDomain::UnpackLM(FEElement &el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N * 6);
	for (int i = 0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3 * i    ] = id[m_dofU[0]];
		lm[3 * i + 1] = id[m_dofU[1]];
		lm[3 * i + 2] = id[m_dofU[2]];

		// rigid rotational dofs
		lm[3 * N + 3 * i    ] = id[m_dofR[0]];
		lm[3 * N + 3 * i + 1] = id[m_dofR[1]];
		lm[3 * N + 3 * i + 2] = id[m_dofR[2]];
	}
}

//-----------------------------------------------------------------------------
void FEDeformableSpringDomain::Activate()
{
	for (int i = 0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.set_active(m_dofU[0]);
				node.set_active(m_dofU[1]);
				node.set_active(m_dofU[2]);
			}
		}
	}

	// calculate the intitial spring length
	m_L0 = InitialLength();
}

//-----------------------------------------------------------------------------
double FEDeformableSpringDomain::InitialLength()
{
	FEMesh& mesh = *m_pMesh;

	double L = 0.0;
	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// get the nodes
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r1 = n1.m_r0;
		vec3d& r2 = n2.m_r0;

		L += (r2 - r1).norm();
	}
	return L;
}

//-----------------------------------------------------------------------------
double FEDeformableSpringDomain::CurrentLength()
{
	FEMesh& mesh = *m_pMesh;

	double L = 0.0;
	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// get the nodes
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r1 = n1.m_rt;
		vec3d& r2 = n2.m_rt;

		L += (r2 - r1).norm();
	}
	return L;
}

//-----------------------------------------------------------------------------
//! Calculates the forces due to discrete elements (i.e. springs)

void FEDeformableSpringDomain::InternalForces(FEGlobalVector& R)
{
	FEMesh& mesh = *m_pMesh;

	vector<double> fe(6);
	vec3d u1, u2;

	vector<int> en(2), lm(6);

	// calculate current length
	double L = CurrentLength();
	double DL = L - m_L0;

	// calculate force
	double F = m_pMat->force(DL);

	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// get the nodes
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r01 = n1.m_r0;
		vec3d& r02 = n2.m_r0;
		vec3d& rt1 = n1.m_rt;
		vec3d& rt2 = n2.m_rt;

		vec3d e = rt2 - rt1; e.unit();

		// calculate spring lengths
		double L0 = (r02 - r01).norm();
		double Lt = (rt2 - rt1).norm();
		double DL = Lt - L0;

		// set up the force vector
		fe[0] = F*e.x;
		fe[1] = F*e.y;
		fe[2] = F*e.z;
		fe[3] = -F*e.x;
		fe[4] = -F*e.y;
		fe[5] = -F*e.z;

		// setup the node vector
		en[0] = el.m_node[0];
		en[1] = el.m_node[1];

		// set up the LM vector
		UnpackLM(el, lm);

		// assemble element
		R.Assemble(en, lm, fe);
	}

	if (m_kbend > 0)
	{
/*		double eps = m_kbend;
		lm.resize(3);
		en.resize(1);
		fe.resize(3);
		int NN = Nodes();
		for (int i = 1; i<NN - 1; ++i)
		{
			int i0 = i - 1;
			int i1 = i + 1;

			vec3d xi = Node(i).m_rt;
			vec3d x0 = Node(i0).m_rt;
			vec3d x1 = Node(i1).m_rt;

			vec3d r = xi - x0;
			vec3d s = x1 - x0; s.unit();
			vec3d d = r - s*(r*s);

			fe[0] = -eps*d.x;
			fe[1] = -eps*d.y;
			fe[2] = -eps*d.z;

			en[0] = m_Node[i];
			lm[0] = Node(i).m_ID[m_dofX];
			lm[1] = Node(i).m_ID[m_dofY];
			lm[2] = Node(i).m_ID[m_dofZ];
			R.Assemble(en, lm, fe);
		}
*/
		double eps = m_kbend;
		lm.resize(3);
		en.resize(1);
		fe.resize(3);
		int NN = Nodes();
		for (int i = 1; i<NN - 1; ++i)
		{
			int i0 = i - 1;
			int i1 = i + 1;

			vec3d xi = Node(i).m_rt;
			vec3d x0 = Node(i0).m_rt;
			vec3d x1 = Node(i1).m_rt;

			vec3d d = xi - (x0 + x1)*0.5;

			fe[0] = -eps*d.x;
			fe[1] = -eps*d.y;
			fe[2] = -eps*d.z;

			en[0] = m_Node[i];
			lm[0] = Node(i).m_ID[m_dofU[0]];
			lm[1] = Node(i).m_ID[m_dofU[1]];
			lm[2] = Node(i).m_ID[m_dofU[2]];
			R.Assemble(en, lm, fe);
		}

	}

	if (m_kstab > 0)
	{
		double eps = m_kstab;
		lm.resize(6);
		en.resize(2);
		fe.resize(6);
		int NE = Elements();
		for (int i=0; i<NE; ++i)
		{
			FEDiscreteElement& el = Element(i);
			en[0] = el.m_node[0];
			en[1] = el.m_node[1];

			// get the nodes
			FENode& n0 = mesh.Node(en[0]);
			FENode& n1 = mesh.Node(en[1]);

			lm[0] = n0.m_ID[m_dofU[0]];
			lm[1] = n0.m_ID[m_dofU[1]];
			lm[2] = n0.m_ID[m_dofU[2]];

			lm[3] = n1.m_ID[m_dofU[0]];
			lm[4] = n1.m_ID[m_dofU[1]];
			lm[5] = n1.m_ID[m_dofU[2]];

			vec3d ei = n1.m_rt - n0.m_rt;

			fe[0] =  eps*ei.x;
			fe[1] =  eps*ei.y;
			fe[2] =  eps*ei.z;
			fe[3] = -eps*ei.x;
			fe[4] = -eps*ei.y;
			fe[5] = -eps*ei.z;

			R.Assemble(en, lm, fe);
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the discrete element stiffness

void FEDeformableSpringDomain::StiffnessMatrix(FELinearSystem& LS)
{
	FEMesh& mesh = *m_pMesh;

	// calculate current length
	double L = CurrentLength();
	double DL = L - m_L0;

	// evaluate the stiffness
	double F = m_pMat->force(DL);
	double E = m_pMat->stiffness(DL);

	FEElementMatrix ke;
	ke.resize(6, 6);
	ke.zero();
	vector<int> en(2), lm(6);

	// loop over all discrete elements
	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// get the nodes of the element
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r01 = n1.m_r0;
		vec3d& r02 = n2.m_r0;
		vec3d& rt1 = n1.m_rt;
		vec3d& rt2 = n2.m_rt;

		vec3d e = rt2 - rt1; e.unit();

		// calculate nodal displacements
		vec3d u1 = rt1 - r01;
		vec3d u2 = rt2 - r02;

		// calculate spring lengths
		double L0 = (r02 - r01).norm();
		double Lt = (rt2 - rt1).norm();
		double DL = Lt - L0;


		if (Lt == 0) { F = 0; Lt = 1; e = vec3d(1, 1, 1); }

		double A[3][3] = { 0 };
		A[0][0] = ((E - F / Lt)*e.x*e.x + F / Lt);
		A[1][1] = ((E - F / Lt)*e.y*e.y + F / Lt);
		A[2][2] = ((E - F / Lt)*e.z*e.z + F / Lt);

		A[0][1] = A[1][0] = (E - F / Lt)*e.x*e.y;
		A[1][2] = A[2][1] = (E - F / Lt)*e.y*e.z;
		A[0][2] = A[2][0] = (E - F / Lt)*e.x*e.z;

		ke[0][0] = A[0][0]; ke[0][1] = A[0][1]; ke[0][2] = A[0][2];
		ke[1][0] = A[1][0]; ke[1][1] = A[1][1]; ke[1][2] = A[1][2];
		ke[2][0] = A[2][0]; ke[2][1] = A[2][1]; ke[2][2] = A[2][2];

		ke[0][3] = -A[0][0]; ke[0][4] = -A[0][1]; ke[0][5] = -A[0][2];
		ke[1][3] = -A[1][0]; ke[1][4] = -A[1][1]; ke[1][5] = -A[1][2];
		ke[2][3] = -A[2][0]; ke[2][4] = -A[2][1]; ke[2][5] = -A[2][2];

		ke[3][0] = -A[0][0]; ke[3][1] = -A[0][1]; ke[3][2] = -A[0][2];
		ke[4][0] = -A[1][0]; ke[4][1] = -A[1][1]; ke[4][2] = -A[1][2];
		ke[5][0] = -A[2][0]; ke[5][1] = -A[2][1]; ke[5][2] = -A[2][2];

		ke[3][3] = A[0][0]; ke[3][4] = A[0][1]; ke[3][5] = A[0][2];
		ke[4][3] = A[1][0]; ke[4][4] = A[1][1]; ke[4][5] = A[1][2];
		ke[5][3] = A[2][0]; ke[5][4] = A[2][1]; ke[5][5] = A[2][2];

		// setup the node vector
		en[0] = el.m_node[0];
		en[1] = el.m_node[1];

		// set up the LM vector
		UnpackLM(el, lm);

		// assemble the element into the global system
		ke.SetNodes(en);
		ke.SetIndices(lm);
		LS.Assemble(ke);
	}

	// Add Bending stiffness
	if (m_kbend > 0)
	{
		double eps = m_kbend;
		vector<int> lmi(3);
		vector<int>	lmj(9);
		en.resize(3);
		int NN = Nodes();
		for (int i = 1; i<NN - 1; ++i)
		{
			int i0 = i - 1;
			int i1 = i + 1;

			ke.resize(3, 9);
			ke.zero();
			ke[0][0] = -eps; ke[0][3] = 0.5*eps; ke[0][6] = 0.5*eps;
			ke[1][1] = -eps; ke[1][4] = 0.5*eps; ke[1][7] = 0.5*eps;
			ke[2][2] = -eps; ke[2][5] = 0.5*eps; ke[2][8] = 0.5*eps;

			vector<int>& IDi = Node(i).m_ID;
			vector<int>& ID0 = Node(i0).m_ID;
			vector<int>& ID1 = Node(i1).m_ID;

			lmi[0] = IDi[m_dofU[0]];
			lmi[1] = IDi[m_dofU[1]];
			lmi[2] = IDi[m_dofU[2]];

			lmj[0] = IDi[m_dofU[0]];
			lmj[1] = IDi[m_dofU[1]];
			lmj[2] = IDi[m_dofU[2]];
			lmj[3] = ID0[m_dofU[0]];
			lmj[4] = ID0[m_dofU[1]];
			lmj[5] = ID0[m_dofU[2]];
			lmj[6] = ID1[m_dofU[0]];
			lmj[7] = ID1[m_dofU[1]];
			lmj[8] = ID1[m_dofU[2]];

			ke.SetNodes(en);
			ke.SetIndices(lmi, lmj);
			LS.Assemble(ke);
		}
	}

	if (m_kstab > 0)
	{
		double eps = m_kstab;
		lm.resize(6);
		en.resize(2);
		ke.resize(6,6); ke.zero();
		ke[0][0] = ke[1][1] = ke[2][2] =  eps;
		ke[3][0] = ke[4][1] = ke[5][2] = -eps;
		ke[0][3] = ke[1][4] = ke[2][5] = -eps;
		ke[3][3] = ke[4][4] = ke[5][5] =  eps;

		int NE = Elements();
		for (int i = 0; i<NE; ++i)
		{
			FEDiscreteElement& el = Element(i);
			en[0] = el.m_node[0];
			en[1] = el.m_node[1];

			// get the nodes
			FENode& n0 = mesh.Node(en[0]);
			FENode& n1 = mesh.Node(en[1]);

			lm[0] = n0.m_ID[m_dofU[0]];
			lm[1] = n0.m_ID[m_dofU[1]];
			lm[2] = n0.m_ID[m_dofU[2]];
			lm[3] = n1.m_ID[m_dofU[0]];
			lm[4] = n1.m_ID[m_dofU[1]];
			lm[5] = n1.m_ID[m_dofU[2]];

			// assemble the element into the global system
			ke.SetNodes(en);
			ke.SetIndices(lm);
			LS.Assemble(ke);
		}
	}
}

//=============================================================================

//-----------------------------------------------------------------------------
FEDeformableSpringDomain2::FEDeformableSpringDomain2(FEModel* pfem) : FEDiscreteDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem), m_dofR(pfem), m_dof(pfem)
{
	m_pMat = 0;

	m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
}

//-----------------------------------------------------------------------------
//! get the total dofs
const FEDofList& FEDeformableSpringDomain2::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEDeformableSpringDomain2::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FESpringMaterial*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
// Only two nodes contribute to this spring
void FEDeformableSpringDomain2::UnpackLM(FEElement &el, vector<int>& lm)
{
	int N = Nodes();
	lm.resize(2 * 6);
	for (int i = 0; i<2; ++i)
	{
		int n = (i==0? 0 : N-1);
		FENode& node = Node(n);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3 * i    ] = id[m_dofU[0]];
		lm[3 * i + 1] = id[m_dofU[1]];
		lm[3 * i + 2] = id[m_dofU[2]];

		// rigid rotational dofs
		lm[3  + 3 * i    ] = id[m_dofR[0]];
		lm[3  + 3 * i + 1] = id[m_dofR[1]];
		lm[3  + 3 * i + 2] = id[m_dofR[2]];
	}
}

//-----------------------------------------------------------------------------
bool FEDeformableSpringDomain2::Init()
{
	if (FEDiscreteDomain::Init() == false) return false;

	// initialize node data
	int NN = Nodes();
	m_nodeData.resize(NN);
	for (int i=0; i<NN; ++i) m_nodeData[i].banchor = false;

	// anchor first and last
	m_nodeData[0].banchor = true;
	m_nodeData[NN-1].banchor = true;

	return true;
}

//-----------------------------------------------------------------------------
void FEDeformableSpringDomain2::Activate()
{
	int N = Nodes();
	for (int i = 0; i<2; ++i)
	{
		int n = (i == 0 ? 0 : N - 1);
		FENode& node = Node(n);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.set_active(m_dofU[0]);
				node.set_active(m_dofU[1]);
				node.set_active(m_dofU[2]);
			}
		}
	}

	// Not sure if this necessary, but let's just do it to be sure
	UpdateNodes();

	// calculate the spring lengths
	m_L0 = InitialLength();
	m_Lt = CurrentLength();
}

//-----------------------------------------------------------------------------
double FEDeformableSpringDomain2::InitialLength()
{
	FEMesh& mesh = *m_pMesh;

	double L = 0.0;
	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// get the nodes
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r1 = n1.m_r0;
		vec3d& r2 = n2.m_r0;

		double DL = (r2 - r1).norm();
		L += DL;
	}
	return L;
}

//-----------------------------------------------------------------------------
double FEDeformableSpringDomain2::CurrentLength()
{
	FEMesh& mesh = *m_pMesh;

	double L = 0.0;
	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		// get the discrete element
		FEDiscreteElement& el = m_Elem[i];

		// get the nodes
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r1 = n1.m_rt;
		vec3d& r2 = n2.m_rt;

		double DL = (r2 - r1).norm();
		L += DL;
	}
	return L;
}

//-----------------------------------------------------------------------------
//! Calculates the forces due to discrete elements (i.e. springs)

void FEDeformableSpringDomain2::InternalForces(FEGlobalVector& R)
{
	FEMesh& mesh = *m_pMesh;

	vector<double> fe(6);
	vec3d u1, u2;

	vector<int> en(2), lm(6);

	// calculate length increment
	double DL = m_Lt - m_L0;

	// calculate force
	double F = m_pMat->force(DL);

	int N = Elements();
	for (size_t i = 0; i<2; ++i)
	{
		int n = (i==0? 0 : N-1);
		double sign = (i == 0 ? 1 : -1);

		// get the discrete element
		FEDiscreteElement& el = Element(n);

		// get the nodes
		FENode& n1 = mesh.Node(el.m_node[0]);
		FENode& n2 = mesh.Node(el.m_node[1]);

		// get the nodal positions
		vec3d& r01 = n1.m_r0;
		vec3d& r02 = n2.m_r0;
		vec3d& rt1 = n1.m_rt;
		vec3d& rt2 = n2.m_rt;

		vec3d e = rt2 - rt1; e.unit();

		// set up the force vector
		fe[i*3    ] = sign*F*e.x;
		fe[i*3 + 1] = sign*F*e.y;
		fe[i*3 + 2] = sign*F*e.z;

		// setup the node vector
		en[i] = el.m_node[i];

		// set up the LM vector
		if (i==0) UnpackLM(el, lm);
	}

	// assemble element
	R.Assemble(en, lm, fe);
}

//-----------------------------------------------------------------------------
//! Calculates the discrete element stiffness

void FEDeformableSpringDomain2::StiffnessMatrix(FELinearSystem& LS)
{
	FEMesh& mesh = *m_pMesh;

	FEElementMatrix ke;
	ke.resize(6, 6);
	ke.zero();

	vector<int> en(2), lm(6);

	// get the nodes of the element
	int NN = Nodes();
	FENode& n1 = Node(0);
	FENode& n2 = Node(NN-1);

	// get the nodal positions
	vec3d& r01 = n1.m_r0;
	vec3d& r02 = n2.m_r0;
	vec3d& rt1 = n1.m_rt;
	vec3d& rt2 = n2.m_rt;

	vec3d e = rt2 - rt1; e.unit();

	// calculate spring lengths
	double L0 = (r02 - r01).norm();
	double Lt = (rt2 - rt1).norm();
	double DL = Lt - L0;

	// evaluate the stiffness
	double F = m_pMat->force(DL);
	double E = m_pMat->stiffness(DL);

	if (Lt == 0) { F = 0; Lt = 1; e = vec3d(1, 1, 1); }

	double A[3][3] = { 0 };
	A[0][0] = ((E - F / Lt)*e.x*e.x + F / Lt);
	A[1][1] = ((E - F / Lt)*e.y*e.y + F / Lt);
	A[2][2] = ((E - F / Lt)*e.z*e.z + F / Lt);

	A[0][1] = A[1][0] = (E - F / Lt)*e.x*e.y;
	A[1][2] = A[2][1] = (E - F / Lt)*e.y*e.z;
	A[0][2] = A[2][0] = (E - F / Lt)*e.x*e.z;

	ke[0][0] = A[0][0]; ke[0][1] = A[0][1]; ke[0][2] = A[0][2];
	ke[1][0] = A[1][0]; ke[1][1] = A[1][1]; ke[1][2] = A[1][2];
	ke[2][0] = A[2][0]; ke[2][1] = A[2][1]; ke[2][2] = A[2][2];

	ke[0][3] = -A[0][0]; ke[0][4] = -A[0][1]; ke[0][5] = -A[0][2];
	ke[1][3] = -A[1][0]; ke[1][4] = -A[1][1]; ke[1][5] = -A[1][2];
	ke[2][3] = -A[2][0]; ke[2][4] = -A[2][1]; ke[2][5] = -A[2][2];

	ke[3][0] = -A[0][0]; ke[3][1] = -A[0][1]; ke[3][2] = -A[0][2];
	ke[4][0] = -A[1][0]; ke[4][1] = -A[1][1]; ke[4][2] = -A[1][2];
	ke[5][0] = -A[2][0]; ke[5][1] = -A[2][1]; ke[5][2] = -A[2][2];

	ke[3][3] = A[0][0]; ke[3][4] = A[0][1]; ke[3][5] = A[0][2];
	ke[4][3] = A[1][0]; ke[4][4] = A[1][1]; ke[4][5] = A[1][2];
	ke[5][3] = A[2][0]; ke[5][4] = A[2][1]; ke[5][5] = A[2][2];

	// setup the node vector
	en[0] = m_Node[0];
	en[1] = m_Node[NN-1];

	// set up the LM vector
	UnpackLM(Element(0), lm);

	// assemble the element into the global system
	ke.SetNodes(en);
	ke.SetIndices(lm);
	LS.Assemble(ke);
}

//-----------------------------------------------------------------------------
void FEDeformableSpringDomain2::Update(const FETimeInfo& tp)
{
	// update wire partition and nodal positions
	UpdateNodes();
}

//-----------------------------------------------------------------------------
void FEDeformableSpringDomain2::SetNodePosition(int node, const vec3d& r)
{
	FENode& nd = Node(node);
	nd.m_rt = r;
	vec3d u = nd.m_rt - nd.m_r0;
	nd.set_vec3d(m_dofU[0], m_dofU[1], m_dofU[2], u);
}

//-----------------------------------------------------------------------------
// This functions distributes the nodes that are not anchored evenly between anchor
// points.
void FEDeformableSpringDomain2::UpdateNodes()
{
	// make sure we have enough nodes
	// We need more than 2 nodes since the outer nodes are always anchored.
	int NN = Nodes();
	if (NN <= 2) return;

	// find wire segments
	int n0 = 0, n1 = 1;
	while (n0 < NN-1)
	{
		if (m_nodeData[n1].banchor)
		{
			vec3d r0 = Node(n0).m_rt;
			vec3d r1 = Node(n1).m_rt;

			for (int n = n0+1; n<=n1-1; ++n)
			{
				assert(m_nodeData[n].banchor == false);

				double w = (double) (n - n0) / (double)(n1 - n0);

				FENode& nd = Node(n);
				nd.m_rt = r0 + (r1 - r0)*w;
				vec3d u = nd.m_rt - nd.m_r0;
				nd.set_vec3d(m_dofU[0], m_dofU[1], m_dofU[2], u);
			}

			n0 = n1;
		}
		n1++;
	}

	// re-calculate current length
	m_Lt = CurrentLength();
}

//-----------------------------------------------------------------------------
//! Anchor (or release) a node
void FEDeformableSpringDomain2::AnchorNode(int node, bool banchor)
{
	m_nodeData[node].banchor = banchor;
}

//-----------------------------------------------------------------------------
vec3d FEDeformableSpringDomain2::NodalForce(int node)
{
	int NN = Nodes();
	if ((node <= 0) || (node >= NN -1)) return vec3d(0,0,0);

	vec3d r = Node(node).m_rt;
	vec3d rm = Node(node-1).m_rt;
	vec3d rp = Node(node+1).m_rt;

	double F = SpringForce();

	vec3d A = rm - r; A.unit(); A *= F;
	vec3d B = rp - r; B.unit(); B *= F;

	vec3d D = A + B;

	return D;
}

//-----------------------------------------------------------------------------
//! get net spring force
double FEDeformableSpringDomain2::SpringForce()
{
	return m_pMat->force(m_Lt - m_L0);
}

//-----------------------------------------------------------------------------
//! tangent
vec3d FEDeformableSpringDomain2::Tangent(int node)
{
	int NN = Nodes();
	if (NN < 2) return vec3d(0,0,0);
	if (node == 0)
	{
		vec3d t = Node(node+1).m_rt - Node(node).m_rt;
		t.unit();
		return t;
	}
	else if (node == NN-1)
	{
		vec3d t = Node(node).m_rt - Node(node - 1).m_rt;
		t.unit();
		return t;
	}
	else
	{
		vec3d t = Node(node + 1).m_rt - Node(node - 1).m_rt;
		t.unit();
		return t;
	}
}
