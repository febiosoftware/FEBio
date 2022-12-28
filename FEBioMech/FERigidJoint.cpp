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
#include "FERigidJoint.h"
#include "FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FEMaterial.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidJoint, FERigidConnector);
	ADD_PARAMETER(m_laugon , "laugon"   );
	ADD_PARAMETER(m_atol   , "tolerance");
	ADD_PARAMETER(m_eps    , "penalty"  );
	ADD_PARAMETER(m_q0     , "joint"    );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidJoint::FERigidJoint(FEModel* pfem) : FERigidConnector(pfem)
{
	static int count = 1;
	m_nID = count++;

	m_laugon = 1;		// Augmented Lagrangian by default for backward compatibility
	m_eps = 0.0;
	m_atol = 0.01;

	m_rbA = m_rbB = nullptr;
}

//-----------------------------------------------------------------------------
//! destructor
FERigidJoint::~FERigidJoint()
{

}

//-----------------------------------------------------------------------------
bool FERigidJoint::Init()
{
	// reset force
	m_F = vec3d(0,0,0);

	// base class first
	if (FERigidConnector::Init() == false) return false;

	// initialize relative joint positions
	m_qa0 = m_q0 - m_rbA->m_r0;
	m_qb0 = m_q0 - m_rbB->m_r0;

	// we make we have a non-zero penalty for penalty and auglag method
	if ((m_laugon != 2) && (m_eps == 0.0)) return false;

	return true;
}

//-----------------------------------------------------------------------------
// allocate equations
int FERigidJoint::InitEquations(int neq)
{
	m_LM.resize(3, -1);
	if (m_laugon == 2)
	{
		// we allocate three equations
		m_LM[0] = neq;
		m_LM[1] = neq + 1;
		m_LM[2] = neq + 2;
		return 3;
	}
	else return 0;
}

//-----------------------------------------------------------------------------
// Build the matrix profile
void FERigidJoint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> lm;
	UnpackLM(lm);

	// add it to the pile
	M.build_add(lm);
}

//-----------------------------------------------------------------------------
void FERigidJoint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	quatd Qa = RBa.GetRotation();
	vec3d qa = Qa*m_qa0;

	quatd Qb = RBb.GetRotation();
	vec3d qb = Qb*m_qb0;

	vec3d ra = RBa.m_rt;
	vec3d rb = RBb.m_rt;
	vec3d c = ra + qa - rb - qb;

	mat3da yaT(qa);
	mat3da ybT(qb);

	vec3d Fa =  m_F;
	vec3d Fb = -m_F;
	vec3d Ma =  yaT*m_F;
	vec3d Mb = -ybT*m_F;

	vector<double> fe(15, 0.0);
	fe[ 0] = -Fa.x;
	fe[ 1] = -Fa.y;
	fe[ 2] = -Fa.z;
	fe[ 3] = -Ma.x;
	fe[ 4] = -Ma.y;
	fe[ 5] = -Ma.z;
	fe[ 6] = -Fb.x;
	fe[ 7] = -Fb.y;
	fe[ 8] = -Fb.z;
	fe[ 9] = -Mb.x;
	fe[10] = -Mb.y;
	fe[11] = -Mb.z;
	fe[12] = -c.x;
	fe[13] = -c.y;
	fe[14] = -c.z;

	vector<int> lm;
	UnpackLM(lm);

	R.Assemble(lm, fe);

	RBa.m_Fr -= vec3d(fe[0], fe[1], fe[2]);
	RBa.m_Mr -= vec3d(fe[3], fe[4], fe[5]);
	RBb.m_Fr -= vec3d(fe[6], fe[7], fe[8]);
	RBb.m_Mr -= vec3d(fe[9], fe[10], fe[11]);
}

//-----------------------------------------------------------------------------
void FERigidJoint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	vec3d a = m_qa0;
	quatd Qa = RBa.GetRotation();
	Qa.RotateVector(a);

	vec3d b = m_qb0;
	quatd Qb = m_rbB->GetRotation();
	Qb.RotateVector(b);

	FEElementMatrix ke;

	if (m_laugon != 2)
	{
		ke.resize(12, 12);
		ke.zero();

		mat3d y1;
		y1[0][0] = 0; y1[0][1] = a.z; y1[0][2] = -a.y;
		y1[1][0] = -a.z; y1[1][1] = 0; y1[1][2] = a.x;
		y1[2][0] = a.y; y1[2][1] = -a.x; y1[2][2] = 0;

		mat3d y2;
		y2[0][0] = 0; y2[0][1] = b.z; y2[0][2] = -b.y;
		y2[1][0] = -b.z; y2[1][1] = 0; y2[1][2] = b.x;
		y2[2][0] = b.y; y2[2][1] = -b.x; y2[2][2] = 0;

		mat3d y11, y12, y22;

		for (int j = 0; j < 3; ++j)
		{
			y11[j][0] = y1[0][j] * y1[0][0] + y1[1][j] * y1[1][0] + y1[2][j] * y1[2][0];
			y11[j][1] = y1[0][j] * y1[0][1] + y1[1][j] * y1[1][1] + y1[2][j] * y1[2][1];
			y11[j][2] = y1[0][j] * y1[0][2] + y1[1][j] * y1[1][2] + y1[2][j] * y1[2][2];

			y12[j][0] = y1[0][j] * y2[0][0] + y1[1][j] * y2[1][0] + y1[2][j] * y2[2][0];
			y12[j][1] = y1[0][j] * y2[0][1] + y1[1][j] * y2[1][1] + y1[2][j] * y2[2][1];
			y12[j][2] = y1[0][j] * y2[0][2] + y1[1][j] * y2[1][2] + y1[2][j] * y2[2][2];

			y22[j][0] = y2[0][j] * y2[0][0] + y2[1][j] * y2[1][0] + y2[2][j] * y2[2][0];
			y22[j][1] = y2[0][j] * y2[0][1] + y2[1][j] * y2[1][1] + y2[2][j] * y2[2][1];
			y22[j][2] = y2[0][j] * y2[0][2] + y2[1][j] * y2[1][2] + y2[2][j] * y2[2][2];
		}

		ke[0][0] = ke[1][1] = ke[2][2] =  1;
		ke[0][6] = ke[1][7] = ke[2][8] = -1;
		ke[6][0] = ke[7][1] = ke[8][2] = -1;
		ke[6][6] = ke[7][7] = ke[8][8] =  1;

		ke[0][3] = y1[0][0]; ke[0][4] = y1[0][1]; ke[0][5] = y1[0][2];
		ke[1][3] = y1[1][0]; ke[1][4] = y1[1][1]; ke[1][5] = y1[1][2];
		ke[2][3] = y1[2][0]; ke[2][4] = y1[2][1]; ke[2][5] = y1[2][2];

		ke[0][9] = -y2[0][0]; ke[0][10] = -y2[0][1]; ke[0][11] = -y2[0][2];
		ke[1][9] = -y2[1][0]; ke[1][10] = -y2[1][1]; ke[1][11] = -y2[1][2];
		ke[2][9] = -y2[2][0]; ke[2][10] = -y2[2][1]; ke[2][11] = -y2[2][2];

		ke[3][0] = y1[0][0]; ke[3][1] = y1[1][0]; ke[3][2] = y1[2][0];
		ke[4][0] = y1[0][1]; ke[4][1] = y1[1][1]; ke[4][2] = y1[2][1];
		ke[5][0] = y1[0][2]; ke[5][1] = y1[1][2]; ke[5][2] = y1[2][2];

		ke[3][3] = y11[0][0]; ke[3][4] = y11[0][1]; ke[3][5] = y11[0][2];
		ke[4][3] = y11[1][0]; ke[4][4] = y11[1][1]; ke[4][5] = y11[1][2];
		ke[5][3] = y11[2][0]; ke[5][4] = y11[2][1]; ke[5][5] = y11[2][2];

		ke[3][6] = -y1[0][0]; ke[3][7] = -y1[1][0]; ke[3][8] = -y1[2][0];
		ke[4][6] = -y1[0][1]; ke[4][7] = -y1[1][1]; ke[4][8] = -y1[2][1];
		ke[5][6] = -y1[0][2]; ke[5][7] = -y1[1][2]; ke[5][8] = -y1[2][2];

		ke[3][9] = -y12[0][0]; ke[3][10] = -y12[0][1]; ke[3][11] = -y12[0][2];
		ke[4][9] = -y12[1][0]; ke[4][10] = -y12[1][1]; ke[4][11] = -y12[1][2];
		ke[5][9] = -y12[2][0]; ke[5][10] = -y12[2][1]; ke[5][11] = -y12[2][2];

		ke[6][3] = -y1[0][0]; ke[6][4] = -y1[0][1]; ke[6][5] = -y1[0][2];
		ke[7][3] = -y1[1][0]; ke[7][4] = -y1[1][1]; ke[7][5] = -y1[1][2];
		ke[8][3] = -y1[2][0]; ke[8][4] = -y1[2][1]; ke[8][5] = -y1[2][2];

		ke[6][9] = y2[0][0]; ke[6][10] = y2[0][1]; ke[6][11] = y2[0][2];
		ke[7][9] = y2[1][0]; ke[7][10] = y2[1][1]; ke[7][11] = y2[1][2];
		ke[8][9] = y2[2][0]; ke[8][10] = y2[2][1]; ke[8][11] = y2[2][2];

		ke[ 9][0] = -y2[0][0]; ke[ 9][1] = -y2[1][0]; ke[ 9][2] = -y2[2][0];
		ke[10][0] = -y2[0][1]; ke[10][1] = -y2[1][1]; ke[10][2] = -y2[2][1];
		ke[11][0] = -y2[0][2]; ke[11][1] = -y2[1][2]; ke[11][2] = -y2[2][2];

		ke[ 9][3] = -y12[0][0]; ke[ 9][4] = -y12[1][0]; ke[ 9][5] = -y12[2][0];
		ke[10][3] = -y12[0][1]; ke[10][4] = -y12[1][1]; ke[10][5] = -y12[2][1];
		ke[11][3] = -y12[0][2]; ke[11][4] = -y12[1][2]; ke[11][5] = -y12[2][2];

		ke[ 9][6] = y2[0][0]; ke[ 9][7] = y2[1][0]; ke [9][8] = y2[2][0];
		ke[10][6] = y2[0][1]; ke[10][7] = y2[1][1]; ke[10][8] = y2[2][1];
		ke[11][6] = y2[0][2]; ke[11][7] = y2[1][2]; ke[11][8] = y2[2][2];

		ke[ 9][9] = y22[0][0]; ke[ 9][10] = y22[0][1]; ke[ 9][11] = y22[0][2];
		ke[10][9] = y22[1][0]; ke[10][10] = y22[1][1]; ke[10][11] = y22[1][2];
		ke[11][9] = y22[2][0]; ke[11][10] = y22[2][1]; ke[11][11] = y22[2][2];

		// scale by penalty factor
		ke *= m_eps;
	}
	else
	{
		ke.resize(15, 15);
		ke.zero();

		mat3dd I(1.0);
		mat3da yaT(a);
		mat3da ybT(b);

		ke.add_symm(0, 12,   I);
		ke.add_symm(3, 12,  yaT);
		ke.add_symm(6, 12,  -I);
		ke.add_symm(9, 12, -ybT);
	}

	// unpack LM
	vector<int> lm;
	UnpackLM(lm);
	ke.SetIndices(lm);

	// assemle into global stiffness matrix
	LS.Assemble(ke);
}

//-----------------------------------------------------------------------------
bool FERigidJoint::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we need to augment
	if (m_laugon == 0) return true;

    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;

	quatd Qa = RBa.GetRotation();
	quatd Qb = RBb.GetRotation();

	vec3d ra = RBa.m_rt;
	vec3d rb = RBb.m_rt;

	vec3d qa = Qa*m_qa0;
	vec3d qb = Qb*m_qb0;

	vec3d c = ra + qa - rb - qb;

	// For Lagrange multipliers we just report the values
	// of the LM and the constraint
	if (m_laugon == 2)
	{
		feLog("\n=== rigid joint # %d:\n", m_nID);
		feLog("\tLagrange m. : %15.7lg, %15.7lg, %15.7lg\n", m_F.x, m_F.y, m_F.z);
		feLog("\tconstraint  : %15.7lg, %15.7lg, %15.7lg\n", c.x, c.y, c.z);
		return true;
	}

	// augmented Lagrangian
	bool bconv = true;

	double normF0 = sqrt(m_L*m_L);

	// calculate trial multiplier
	vec3d Lm = m_L + c*m_eps;

	double normF1 = sqrt(Lm*Lm);

	// check convergence of constraints
	feLog(" rigid joint # %d\n", m_nID);
	feLog("                  CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normF1) > 1e-10) pctn = fabs((normF1 - normF0)/normF1);
	feLog("    force : %15le %15le\n", pctn, m_atol);
	feLog("    gap   : %15le       ***\n", c.norm());
		
	if (pctn >= m_atol)
	{
		bconv = false;

		// update multiplier
		m_L = m_L + c*m_eps;
	
		// update force
		m_F = m_L + c*m_eps;
	}
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FERigidJoint::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
	ar & m_nID;
	ar & m_q0 & m_qa0 & m_qb0;
	ar & m_F & m_L & m_eps & m_atol & m_laugon;
}

//-----------------------------------------------------------------------------
void FERigidJoint::Update()
{
	if (m_laugon != 2)
	{
		FERigidBody& RBa = *m_rbA;
		FERigidBody& RBb = *m_rbB;

		vec3d ra = RBa.m_rt;
		vec3d rb = RBb.m_rt;

		vec3d qa = m_qa0;
		RBa.GetRotation().RotateVector(qa);

		vec3d qb = m_qb0;
		RBb.GetRotation().RotateVector(qb);

		vec3d c = ra + qa - rb - qb;
		m_F = m_L + c*m_eps;
	}
}

//-----------------------------------------------------------------------------
void FERigidJoint::Reset()
{
	m_F = vec3d(0,0,0);
	m_L = vec3d(0,0,0);

	m_qa0 = m_q0 - m_rbA->m_r0;
	m_qb0 = m_q0 - m_rbB->m_r0;
}

//-----------------------------------------------------------------------------
void FERigidJoint::UnpackLM(vector<int>& lm)
{
	// add the dofs of rigid body A
	lm.reserve(15);
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
		lm.push_back(m_LM[0]);
		lm.push_back(m_LM[1]);
		lm.push_back(m_LM[2]);
	}
}

//-----------------------------------------------------------------------------
void FERigidJoint::Update(const std::vector<double>& ui)
{
	if (m_laugon == 2)
	{
		m_F.x += ui[m_LM[0]];
		m_F.y += ui[m_LM[1]];
		m_F.z += ui[m_LM[2]];
	}
}
