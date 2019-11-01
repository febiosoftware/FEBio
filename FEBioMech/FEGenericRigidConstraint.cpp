/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FEGenericRigidConstraint.h"
#include "FERigidBody.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEGenericRigidConstraint, FERigidConnector);
	ADD_PARAMETER(m_nRBa   , "body_a" );
	ADD_PARAMETER(m_nRBb   , "body_b" );
	ADD_PARAMETER(m_q0     , "joint"  );
	ADD_PARAMETER(m_bc[0]  , "fix_x"  );
	ADD_PARAMETER(m_bc[1]  , "fix_y"  );
	ADD_PARAMETER(m_bc[2]  , "fix_z"  );
	ADD_PARAMETER(m_bc[3]  , "fix_Rx" );
	ADD_PARAMETER(m_bc[4]  , "fix_Ry" );
	ADD_PARAMETER(m_bc[5]  , "fix_Rz" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEGenericRigidConstraint::FEGenericRigidConstraint(FEModel* pfem) : FERigidConnector(pfem)
{
	static int count = 1;
	m_nID = count++;

	m_rbA = m_rbB = nullptr;

	m_bc[0] = false;
	m_bc[1] = false;
	m_bc[2] = false;
	m_bc[3] = false;
	m_bc[4] = false;
	m_bc[5] = false;

	m_v0[0] = vec3d(1, 0, 0);
	m_v0[1] = vec3d(0, 1, 0);
	m_v0[2] = vec3d(0, 0, 1);

	for (int i = 0; i < 6; ++i)
	{
		m_EQ[i] = -1;
		m_Lm[i] = 0.0;
		m_Lmp[i] = 0.0;
	}
}

//-----------------------------------------------------------------------------
bool FEGenericRigidConstraint::Init()
{
	// base class first
	if (FERigidConnector::Init() == false) return false;

	// initialize relative joint positions
	m_qa0 = m_q0 - m_rbA->m_r0;
	m_qb0 = m_q0 - m_rbB->m_r0;

	return true;
}

//-----------------------------------------------------------------------------
// allocate equations
int FEGenericRigidConstraint::InitEquations(int neq)
{
	int n = neq;
	for (int i=0; i<6; ++i)
	{
		if (m_bc[i]) m_EQ[i] = n++; else m_EQ[i] = -1;
	}
	return n - neq;
}

//-----------------------------------------------------------------------------
// Build the matrix profile
void FEGenericRigidConstraint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> lm;
	UnpackLM(lm);

	// add it to the pile
	M.build_add(lm);
}

//-----------------------------------------------------------------------------
void FEGenericRigidConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	quatd Qa = RBa.GetRotation();
	vec3d qa = Qa*m_qa0;

	quatd Qb = RBb.GetRotation();
	vec3d qb = Qb*m_qb0;

	vec3d ra = RBa.m_rt;
	vec3d rb = RBb.m_rt;
	vec3d g = ra + qa - rb - qb;

	mat3da ya(-qa); mat3da yaT = -ya;
	mat3da yb(-qb); mat3da ybT = -yb;

	// add translational constraints
	vector<double> fe(18, 0.0);
	for (int i=0; i<3; ++i)
	{ 
		if (m_bc[i])
		{
			double L = m_Lm[i];

			vec3d v = Qa*m_v0[i];
			mat3da VT(v);

			vec3d F = v*L;

			vec3d Fa = F;
			vec3d Fb = -F;
			vec3d Ma =  yaT*F + VT*(g*L);
			vec3d Mb = -ybT*F;

			// body A
			fe[0] -= Fa.x;
			fe[1] -= Fa.y;
			fe[2] -= Fa.z;
			fe[3] -= Ma.x;
			fe[4] -= Ma.y;
			fe[5] -= Ma.z;

			// body B
			fe[6] -= Fb.x;
			fe[7] -= Fb.y;
			fe[8] -= Fb.z;
			fe[9] -= Mb.x;
			fe[10] -= Mb.y;
			fe[11] -= Mb.z;

			// constraint
			fe[12 + i] -= v*g;
		}
	}

	// add rotational constraints
	int LUT[3][2] = { { 1,2 },{ 2, 0 },{ 0, 1 } };
	for (int i = 3; i < 6; ++i)
	{
		if (m_bc[i])
		{
			int K = i - 3;
			int I = LUT[K][0];
			int J = LUT[K][1];
			vec3d va = Qa*m_v0[I];
			vec3d vb = Qb*m_v0[J];

			double L = m_Lm[i];

			mat3da VaT(va);
			mat3da VbT(vb);

			vec3d Ma = VaT*vb*L;
			vec3d Mb = VbT*va*L;

			// body A
			fe[ 3] -= Ma.x;
			fe[ 4] -= Ma.y;
			fe[ 5] -= Ma.z;

			// body B
			fe[ 9] -= Mb.x;
			fe[10] -= Mb.y;
			fe[11] -= Mb.z;

			// constraint
			fe[12 + i] -= va*vb;
		}
	}

	vector<int> lm;
	UnpackLM(lm);

	R.Assemble(lm, fe);

	RBa.m_Fr -= vec3d(fe[0], fe[1], fe[2]);
	RBa.m_Mr -= vec3d(fe[3], fe[4], fe[5]);
	RBb.m_Fr -= vec3d(fe[6], fe[7], fe[8]);
	RBb.m_Mr -= vec3d(fe[9], fe[10], fe[11]);
}

//-----------------------------------------------------------------------------
void FEGenericRigidConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	quatd Qa = RBa.GetRotation();
	vec3d qa = Qa*m_qa0;

	quatd Qb = m_rbB->GetRotation();
	vec3d qb = Qb*m_qb0;

	vec3d ra = RBa.m_rt;
	vec3d rb = RBb.m_rt;
	vec3d g = ra + qa - rb - qb;

	mat3da ya(-qa);
	mat3da yaT = -ya;
	mat3da yb(-qb); 
	mat3da ybT = -yb;

	FEElementMatrix ke;
	ke.resize(18, 18);
	ke.zero();

	// translational constraints
	for (int i = 0; i < 3; ++i)
	{
		if (m_bc[i])
		{
			double L = m_Lm[i];
			vec3d v = Qa*m_v0[i];
			mat3da V(-v);
			mat3da VT = -V;

			ke.add_symm(0, 3, V*L);
			ke.add_symm(0, 12 + i, v);

			ke.add(3, 3, (yaT*V + VT*ya)*L);
			ke.add_symm(3, 6, -VT*L);
			ke.add_symm(3, 9, -(VT*yb)*L);
			ke.add_symm(3, 12 + i, yaT*v + VT*g);

			ke.add_symm(6, 12 + i, -v);

			ke.add_symm(9, 12 + i, -ybT*v);
		}
	}

	// rotational constraints
	int LUT[3][2] = { {1,2}, {2, 0}, {0, 1} };
	for (int i = 3; i < 6; ++i)
	{
		if (m_bc[i])
		{
			int K = i - 3;
			int I = LUT[K][0];
			int J = LUT[K][1];
			vec3d va = Qa*m_v0[I];
			vec3d vb = Qb*m_v0[J];

			mat3da VaT(va);
			mat3da VbT(vb);

			ke.add_symm(3, 12 + i, VaT*vb);
			ke.add_symm(9, 12 + i, VbT*va);
		}
	}

	// unpack LM
	vector<int> lm;
	UnpackLM(lm);
	ke.SetIndices(lm);

	// assemle into global stiffness matrix
	LS.Assemble(ke);
}

//-----------------------------------------------------------------------------
bool FEGenericRigidConstraint::Augment(int naug, const FETimeInfo& tp)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;

	quatd Qa = RBa.GetRotation();
	quatd Qb = RBb.GetRotation();

	vec3d ra = RBa.m_rt;
	vec3d rb = RBb.m_rt;

	vec3d qa = Qa*m_qa0;
	vec3d qb = Qb*m_qb0;

	vec3d g = ra + qa - rb - qb;

	double c[6] = { 0 };

	// translational constraints
	for (int i = 0; i < 3; ++i)
	{
		vec3d v = Qa*m_v0[i];
		c[i] = v*g;
	}

	// rotational constraints
	for (int i = 3; i < 6; ++i)
	{
		int K = i - 3;
		int I = (K + 1) % 3;
		int J = (K + 2) % 3;
		vec3d va = Qa*m_v0[I];
		vec3d vb = Qb*m_v0[J];

		c[i] = va*vb;
	}

	// For Lagrange multipliers we just report the values
	// of the LM and the constraint
	feLog("\n=== rigid connector # %d:\n", m_nID);
	feLog(" Lagrange m. : %13.7lg,%13.7lg,%13.7lg\n", m_Lm[0], m_Lm[1], m_Lm[2]);
	feLog("               %13.7lg,%13.7lg,%13.7lg\n", m_Lm[3], m_Lm[4], m_Lm[5]);
	feLog(" constraint  : %13.7lg,%13.7lg,%13.7lg\n", c[0], c[1], c[2]);
	feLog("               %13.7lg,%13.7lg,%13.7lg\n", c[3], c[4], c[5]);
	return true;
}

//-----------------------------------------------------------------------------
void FEGenericRigidConstraint::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
	ar & m_nID;
	ar & m_q0 & m_qa0 & m_qb0;
	ar & m_Lm & m_Lmp;
}

//-----------------------------------------------------------------------------
void FEGenericRigidConstraint::Update()
{

}

//-----------------------------------------------------------------------------
void FEGenericRigidConstraint::Reset()
{
	assert(false);
}

//-----------------------------------------------------------------------------
void FEGenericRigidConstraint::UnpackLM(vector<int>& lm)
{
	// add the dofs of rigid body A
	lm.reserve(18);
	lm.push_back(m_rbA->m_LM[0]);
	lm.push_back(m_rbA->m_LM[1]);
	lm.push_back(m_rbA->m_LM[2]);
	lm.push_back(m_rbA->m_LM[3]);
	lm.push_back(m_rbA->m_LM[4]);
	lm.push_back(m_rbA->m_LM[5]);

	// add the dofs of rigid body B
	lm.push_back(m_rbB->m_LM[0]);
	lm.push_back(m_rbB->m_LM[1]);
	lm.push_back(m_rbB->m_LM[2]);
	lm.push_back(m_rbB->m_LM[3]);
	lm.push_back(m_rbB->m_LM[4]);
	lm.push_back(m_rbB->m_LM[5]);

	// add the LM equations
	for (int i=0; i<6; ++i) lm.push_back(m_EQ[i]);
}

//-----------------------------------------------------------------------------
void FEGenericRigidConstraint::Update(const std::vector<double>& Ui, const std::vector<double>& ui)
{
	for (int i=0; i<6; ++i)
	{
		if (m_EQ[i] != -1) m_Lm[i] = m_Lmp[i] + Ui[m_EQ[i]] + ui[m_EQ[i]];
	}
}

void FEGenericRigidConstraint::PrepStep()
{
	for (int i = 0; i < 6; ++i)
	{
		m_Lmp[i] = m_Lm[i];
	}
}

void FEGenericRigidConstraint::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui)
{
	for (int i = 0; i < 6; ++i)
	{
		if (m_EQ[i] != -1) Ui[m_EQ[i]] += ui[m_EQ[i]];
	}
}
