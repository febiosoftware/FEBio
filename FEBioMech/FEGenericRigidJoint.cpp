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
#include "FEGenericRigidJoint.h"
#include "FERigidBody.h"
#include <FECore/log.h>
#include <FECore/FEMaterial.h>
#include <FECore/FELinearSystem.h>
#include <FECore/fecore_debug.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEGenericRigidJoint, FERigidConnector);
	ADD_PARAMETER(m_laugon , "laugon" );
	ADD_PARAMETER(m_eps    , "penalty");
	ADD_PARAMETER(m_tol    , "tolerance");
	ADD_PARAMETER(m_q0     , "joint"  );
	ADD_PARAMETER(m_dc[0]  , "prescribe_x" , m_bc[0]);
	ADD_PARAMETER(m_dc[1]  , "prescribe_y" , m_bc[1]);
	ADD_PARAMETER(m_dc[2]  , "prescribe_z" , m_bc[2]);
	ADD_PARAMETER(m_dc[3]  , "prescribe_Rx", m_bc[3]);
	ADD_PARAMETER(m_dc[4]  , "prescribe_Ry", m_bc[4]);
	ADD_PARAMETER(m_dc[5]  , "prescribe_Rz", m_bc[5]);
	ADD_PARAMETER(m_bsymm  , "symmetric_stiffness");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEGenericRigidJoint::FEGenericRigidJoint(FEModel* pfem) : FERigidConnector(pfem)
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

	m_dc[0] = 0.0;
	m_dc[1] = 0.0;
	m_dc[2] = 0.0;
	m_dc[3] = 0.0;
	m_dc[4] = 0.0;
	m_dc[5] = 0.0;

	m_bsymm = false;

	m_laugon = 0;	// default to penalty method
	m_eps = 1.0;
	m_tol = 0.0;

	vec3d a(1, 0, 0);
	vec3d d(0, 1, 0);
	vec3d e1 = a.normalized();
	vec3d e3 = a ^ d; e3.unit();
	vec3d e2 = e3 ^ e1;

	m_v0[0] = e1;
	m_v0[1] = e2;
	m_v0[2] = e3;

	for (int i = 0; i < 6; ++i)
	{
		m_EQ[i] = -1;
		m_Lm[i] = 0.0;
		m_Lmp[i] = 0.0;
	}
}

//-----------------------------------------------------------------------------
bool FEGenericRigidJoint::Init()
{
	// base class first
	if (FERigidConnector::Init() == false) return false;

	// initialize relative joint positions
	m_qa0 = m_q0 - m_rbA->m_r0;
	m_qb0 = m_q0 - m_rbB->m_r0;

	return true;
}

//-----------------------------------------------------------------------------
//! initial position 
vec3d FEGenericRigidJoint::InitialPosition() const
{
	return m_q0;
}

//-----------------------------------------------------------------------------
//! current position
vec3d FEGenericRigidJoint::Position() const
{
	FERigidBody& RBa = *m_rbA;
	quatd Qa = RBa.GetRotation();
	vec3d qa = Qa * m_qa0;
	vec3d ra = RBa.m_rt;
	
	return ra + qa;
}

//-----------------------------------------------------------------------------
// allocate equations
int FEGenericRigidJoint::InitEquations(int neq)
{
	if (m_laugon == 2)
	{
		int n = neq;
		for (int i = 0; i < 6; ++i)
		{
			if (m_bc[i]) m_EQ[i] = n++; else m_EQ[i] = -1;
		}
		return n - neq;
	}
	else return 0;
}

//-----------------------------------------------------------------------------
// Build the matrix profile
void FEGenericRigidJoint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> lm;
	UnpackLM(lm);

	// add it to the pile
	M.build_add(lm);
}

//-----------------------------------------------------------------------------
void FEGenericRigidJoint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
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

	int ndof = (m_laugon == 2 ? 18 : 12);

	// add translational constraints
	vector<double> fe(ndof, 0.0);
	for (int i = 0; i < 3; ++i)
	{
		if (m_bc[i])
		{
			vec3d v = Qa*m_v0[i];
			mat3da VT(v);

			double c = v*g + m_dc[i];

			double L = m_Lm[i];
			if (m_laugon != 2) L += m_eps*c;

			vec3d F = v*L;

			vec3d Fa = F;
			vec3d Fb = -F;
			vec3d Ma = yaT*F + VT*(g*L);
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

			// add constraint for Lagrange multiplier
			if (m_laugon == 2)
			{
				fe[12 + i] -= c;
			}
		}
	}

	// add rotational constraints
	for (int i = 3; i < 6; ++i)
	{
		if (m_bc[i])
		{
			int K = i - 3;
			int I = (K + 1) % 3;
			int J = (K + 2) % 3;
			vec3d va = Qa*m_v0[I];
			vec3d vb = Qb*m_v0[J];

			double c = va*vb - cos(0.5*PI + m_dc[i]);

			double L = m_Lm[i];
			if (m_laugon != 2) L += m_eps*c;

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

			// add constraint for Lagrange multiplier
			if (m_laugon == 2)
			{
				fe[12 + i] -= c;
			}
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
void FEGenericRigidJoint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	if (m_laugon == 2) StiffnessMatrixLM(LS, tp);
	else StiffnessMatrixAL(LS, tp);
}

//-----------------------------------------------------------------------------
void FEGenericRigidJoint::StiffnessMatrixLM(FELinearSystem& LS, const FETimeInfo& tp)
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

	mat3da G(-g);
	mat3da GT = -G;

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
			mat3d V = mat3da(-v);
			mat3d VT = -V;

			ke.add_symm(0, 3, V*L);
			ke.add_symm(0, 12 + i, v);

			if (m_bsymm) ke.add(3, 3, (yaT*V - GT*V).sym()*L);
			else ke.add(3, 3, (yaT*V - GT*V)*L);

			ke.add_symm(3, 6, -VT*L);
			ke.add_symm(3, 9, -(VT*yb)*L);
			ke.add_symm(3, 12 + i, yaT*v + VT*g);

			ke.add_symm(6, 12 + i, -v);

			if (m_bsymm) ke.add(9, 9     , (VT*yb).sym()*L);
			else ke.add(9, 9, (VT*yb)*L);

			ke.add_symm(9, 12 + i, -ybT*v);
		}
	}

	// rotational constraints
	for (int i = 3; i < 6; ++i)
	{
		if (m_bc[i])
		{
			double L = m_Lm[i];

			int K = i - 3;
			int I = (K + 1) % 3;
			int J = (K + 2) % 3;
			vec3d va = Qa*m_v0[I];
			vec3d vb = Qb*m_v0[J];

			mat3da Va(-va), VaT(va);
			mat3da Vb(-vb), VbT(vb);

			if (m_bsymm) ke.add(3, 3, -(VbT*Va).sym()*L);
			else ke.add(3, 3, -(VbT*Va)*L);

			ke.add_symm(3, 9, VaT*Vb*L);
			ke.add_symm(3, 12 + i, VaT*vb);

			if (m_bsymm) ke.add(9, 9, -(VaT*Vb).sym()*L);
			else ke.add(9, 9, -(VaT*Vb)*L);

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
void FEGenericRigidJoint::StiffnessMatrixAL(FELinearSystem& LS, const FETimeInfo& tp)
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

	mat3da G(-g);
	mat3da GT = -G;

	FEElementMatrix ke;
	ke.resize(12, 12);
	ke.zero();

	// translational constraints
	for (int i = 0; i < 3; ++i)
	{
		if (m_bc[i])
		{
			vec3d w = Qa*m_v0[i];
			mat3d W = mat3da(-w);
			mat3d WT = -W;

			double c = w*g + m_dc[i];
			double lam = m_Lm[i] + m_eps*c;

			mat3dd I(1.0);
			matrix Q(12, 3);
			Q.set(0, 0,    I);
			Q.set(3, 0,  yaT);
			Q.set(6, 0,   -I);
			Q.set(9, 0, -ybT);
			matrix Qw = Q*w;

			matrix V(12, 3); V.zero();
			V.set(3, 0, WT);
			matrix Vg = V*g;

			matrix C = Vg + Qw;

			matrix CCT = C*C.transpose();

			matrix QVT = Q*V.transpose();
			matrix VQT = V*Q.transpose();

			ke += CCT*m_eps + (QVT + VQT)*lam;

			if (m_bsymm == false)
			{
				ke.add(3, 3, (GT*W + WT*ya)*(-lam));
				ke.add(9, 9, (WT*yb)*(lam));
			}
		}
	}

	// rotational constraints
	for (int i = 3; i < 6; ++i)
	{
		if (m_bc[i])
		{
			int K = i - 3;
			int I = (K + 1) % 3;
			int J = (K + 2) % 3;
			vec3d va = Qa*m_v0[I];
			vec3d vb = Qb*m_v0[J];

			double c = va*vb - cos(0.5*PI + m_dc[i]);

			double lam = m_Lm[i] + m_eps*c;

			mat3da Va(-va), VaT(va);
			mat3da Vb(-vb), VbT(vb);

			mat3dd Id(1.0);
			matrix N(12, 3); N.zero();
			N.set(3, 0,  Id);
			N.set(9, 0, -Id);
		
			vec3d vv = VaT*vb;
			matrix C = N*vv;

			matrix CCT = C*C.transpose();
			ke += CCT*m_eps;

			if (m_bsymm == false)
			{
				mat3d VV = VaT*Vb;
				matrix WT(3, 12); WT.zero();
				WT.add(0, 3, -VV.transpose());
				WT.add(0, 9, VV);

				matrix JWT = N*WT;
				ke += JWT*lam;
			}
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
bool FEGenericRigidJoint::Augment(int naug, const FETimeInfo& tp)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;

	quatd Qa = RBa.GetRotation();
	quatd Qb = RBb.GetRotation();

	vec3d ra = RBa.m_rt;
	vec3d rb = RBb.m_rt;

	vec3d qa = Qa*m_qa0;
	vec3d qb = Qb*m_qb0;

	// NOTE: This is NOT the gap that was used in the last residual
	//       evaluation. Is that a problem?
	vec3d pa = ra + qa;
	vec3d pb = rb + qb;
	vec3d g = (ra - rb) + (qa - qb);

	double c[6] = { 0 };

	// translational constraints
	for (int i = 0; i < 3; ++i)
	{
		if (m_bc[i])
		{
			vec3d v = Qa*m_v0[i];
			c[i] = v*g + m_dc[i];
		}
	}

	// rotational constraints
	for (int i = 3; i < 6; ++i)
	{
		if (m_bc[i])
		{
			int K = i - 3;
			int I = (K + 1) % 3;
			int J = (K + 2) % 3;
			vec3d va = Qa*m_v0[I];
			vec3d vb = Qb*m_v0[J];

			c[i] = va*vb - cos(0.5*PI + m_dc[i]);
		}
	}

	feLog("\n=== rigid connector # %d:\n", m_nID);

	bool bconv = true;
	double e[6] = { 0 };
	double Lm[6] = { 0 };
	if (m_laugon == 0)
	{
		for (int i = 0; i < 6; ++i)
		{
			Lm[i] = m_eps*c[i];
		}
	}
	else if (m_laugon == 1)
	{
//		if (m_tol > 0.0)
		{
			for (int i = 0; i < 6; ++i)
			{
				if (m_bc[i])
				{
					Lm[i] = m_Lm[i] + m_eps*c[i];

					double F0 = fabs(m_Lm[i]);
					double F1 = fabs(Lm[i]);

					double dl = fabs((F1 - F0) / F1);
					if (m_tol && (dl > m_tol)) bconv = false;

					e[i] = dl;
				}
			}

//			if (bconv == false)
			{
				for (int i = 0; i < 6; ++i)
					m_Lm[i] += m_eps*c[i];
			}
		}
	}
	else
	{
		for (int i = 0; i < 6; ++i) Lm[i] = m_Lm[i];
	}

	feLog(" Lagrange m. : %13.7lg,%13.7lg,%13.7lg\n", Lm[0], Lm[1], Lm[2]);
	feLog("               %13.7lg,%13.7lg,%13.7lg\n", Lm[3], Lm[4], Lm[5]);
	feLog(" constraint  : %13.7lg,%13.7lg,%13.7lg\n", c[0], c[1], c[2]);
	feLog("               %13.7lg,%13.7lg,%13.7lg\n", c[3], c[4], c[5]);
	if (m_laugon == 1)
	{
		feLog(" error       : %13.7lg,%13.7lg,%13.7lg\n", e[0], e[1], e[2]);
		feLog("               %13.7lg,%13.7lg,%13.7lg\n", e[3], e[4], e[5]);
	}

	return bconv;
}

//-----------------------------------------------------------------------------
void FEGenericRigidJoint::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
	ar & m_nID;
	ar & m_q0 & m_qa0 & m_qb0;
	ar & m_Lm & m_Lmp;
}

//-----------------------------------------------------------------------------
void FEGenericRigidJoint::Update()
{

}

//-----------------------------------------------------------------------------
void FEGenericRigidJoint::Reset()
{
	assert(false);
}

//-----------------------------------------------------------------------------
void FEGenericRigidJoint::UnpackLM(vector<int>& lm)
{
	int ndof = (m_laugon == 2 ? 18 : 12);

	// add the dofs of rigid body A
	lm.reserve(ndof);
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
	if (m_laugon == 2)
	{
		for (int i = 0; i < 6; ++i) lm.push_back(m_EQ[i]);
	}
}

//-----------------------------------------------------------------------------
void FEGenericRigidJoint::Update(const std::vector<double>& Ui, const std::vector<double>& ui)
{
	if (m_laugon == 2)
	{
		for (int i = 0; i < 6; ++i)
		{
			if (m_EQ[i] != -1) m_Lm[i] = m_Lmp[i] + Ui[m_EQ[i]] + ui[m_EQ[i]];
		}
	}
}

void FEGenericRigidJoint::PrepStep()
{
	if (m_laugon == 2)
	{
		for (int i = 0; i < 6; ++i)
		{
			m_Lmp[i] = m_Lm[i];
		}
	}
}

void FEGenericRigidJoint::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui)
{
	if (m_laugon == 2)
	{
		for (int i = 0; i < 6; ++i)
		{
			if (m_EQ[i] != -1) Ui[m_EQ[i]] += ui[m_EQ[i]];
		}
	}
}
