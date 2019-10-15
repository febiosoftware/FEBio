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
#include "FELMRigidJoint.h"
#include "FEMechModel.h"
#include "FERigidMaterial.h"
#include "FERigidSystem.h"
#include "FERigidBody.h"
#include <FECore/FELinearSystem.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FELMRigidJoint, FELMConstraint)
	ADD_PARAMETER(m_nRBa   , "body_a"   );
	ADD_PARAMETER(m_nRBb   , "body_b"   );
	ADD_PARAMETER(m_q0     , "joint");
END_FECORE_CLASS();

FELMRigidJoint::FELMRigidJoint(FEModel* fem) : FELMConstraint(fem)
{
	m_nRBa = -1;
	m_nRBb = -1;
	m_q0 = vec3d(0, 0, 0);
}

bool FELMRigidJoint::Init()
{
	// do base class first
	if (FENLConstraint::Init() == false) return false;

	// When the rigid spring is read in, the ID's correspond to the rigid materials.
	// Now we want to make the ID's refer to the rigid body ID's
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_nRBa - 1));
	if (pm == nullptr)
	{
		feLogError("Rigid joint does not connect two rigid bodies\n");
		return false;
	}
	m_nRBa = pm->GetRigidBodyID();

	pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_nRBb - 1));
	if (pm == nullptr)
	{
		feLogError("Rigid connector does not connect two rigid bodies\n");
		return false;
	}
	m_nRBb = pm->GetRigidBodyID();

	// get the actual rigid bodies
	FERigidSystem& rs = *fem.GetRigidSystem();
	m_rbA = rs.Object(m_nRBa);
	m_rbB = rs.Object(m_nRBb);

	// initialize relative joint positions
	m_qa0 = m_q0 - m_rbA->m_r0;
	m_qb0 = m_q0 - m_rbB->m_r0;

	return true;
}

// allocate equations
int FELMRigidJoint::InitEquations(int neq)
{
	// we allocate three equations
	m_LM.resize(3);
	m_LM[0] = neq;
	m_LM[1] = neq + 1;
	m_LM[2] = neq + 2;
	return 3;
}

// The LoadVector function evaluates the "forces" that contribute to the residual of the system
void FELMRigidJoint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	vec3d ra = m_rbA->m_rt;
	vec3d rb = m_rbB->m_rt;

	vec3d a = m_qa0;
	quatd Qa = m_rbA->GetRotation();
	Qa.RotateVector(a);

	vec3d b = m_qb0;
	quatd Qb = m_rbB->GetRotation();
	Qb.RotateVector(b);

	vec3d c = ra + a - rb - b;


	vector<double> fe(15, 0.0);
	fe[ 0] = -m_Lm.x;
	fe[ 1] = -m_Lm.y;
	fe[ 2] = -m_Lm.z;
	fe[ 3] = -a.y*m_Lm.z + a.z*m_Lm.y;
	fe[ 4] = -a.z*m_Lm.x + a.x*m_Lm.z;
	fe[ 5] = -a.x*m_Lm.y + a.y*m_Lm.x;
	fe[ 6] = m_Lm.x;
	fe[ 7] = m_Lm.y;
	fe[ 8] = m_Lm.z;
	fe[ 9] = b.y*m_Lm.z - b.z*m_Lm.y;
	fe[10] = b.z*m_Lm.x - b.x*m_Lm.z;
	fe[11] = b.x*m_Lm.y - b.y*m_Lm.x;
	fe[12] = -c.x;
	fe[13] = -c.y;
	fe[14] = -c.z;

	vector<int> lm;
	UnpackLM(lm);

	R.Assemble(lm, fe);
}

// Evaluates the contriubtion to the stiffness matrix
void FELMRigidJoint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	vec3d a = m_qa0;
	quatd Qa = m_rbA->GetRotation();
	Qa.RotateVector(a);

	mat3d y1;
	y1[0][0] =    0; y1[0][1] =  a.z; y1[0][2] = -a.y;
	y1[1][0] = -a.z; y1[1][1] =    0; y1[1][2] =  a.x;
	y1[2][0] =  a.y; y1[2][1] = -a.x; y1[2][2] =    0;

	vec3d b = m_qb0;
	quatd Qb = m_rbB->GetRotation();
	Qb.RotateVector(b);

	mat3d y2;
	y2[0][0] =    0; y2[0][1] =  b.z; y2[0][2] = -b.y;
	y2[1][0] = -b.z; y2[1][1] =    0; y2[1][2] =  b.x;
	y2[2][0] =  b.y; y2[2][1] = -b.x; y2[2][2] =    0;

	FEElementMatrix ke(15, 15);
	ke.zero();

	ke[0][12] = ke[12][0] = 1;
	ke[1][13] = ke[13][1] = 1;
	ke[2][14] = ke[14][2] = 1;

	ke[3][12] = -y1[0][0]; ke[3][13] = -y1[0][1]; ke[3][14] = -y1[0][2];
	ke[4][12] = -y1[1][0]; ke[4][13] = -y1[1][1]; ke[4][14] = -y1[1][2];
	ke[5][12] = -y1[2][0]; ke[5][13] = -y1[2][1]; ke[5][14] = -y1[2][2];

	ke[6][12] = ke[12][6] = -1;
	ke[7][13] = ke[13][7] = -1;
	ke[8][14] = ke[14][8] = -1;

	ke[ 9][12] = y2[0][0]; ke[ 9][13] = y2[0][1]; ke[ 9][14] = y2[0][2];
	ke[10][12] = y2[1][0]; ke[10][13] = y2[1][1]; ke[10][14] = y2[1][2];
	ke[11][12] = y2[2][0]; ke[11][13] = y2[2][1]; ke[11][14] = y2[2][2];

	ke[12][3] = y1[0][0]; ke[12][4] = y1[0][1]; ke[12][5] = y1[0][2];
	ke[13][3] = y1[1][0]; ke[13][4] = y1[1][1]; ke[13][5] = y1[1][2];
	ke[14][3] = y1[2][0]; ke[14][4] = y1[2][1]; ke[14][5] = y1[2][2];

	ke[12][9] = -y2[0][0]; ke[12][10] = -y2[0][1]; ke[12][11] = -y2[0][2];
	ke[13][9] = -y2[1][0]; ke[13][10] = -y2[1][1]; ke[13][11] = -y2[1][2];
	ke[14][9] = -y2[2][0]; ke[14][10] = -y2[2][1]; ke[14][11] = -y2[2][2];

	vector<int> lm;
	UnpackLM(lm);
	ke.SetIndices(lm);
	LS.Assemble(ke);
}

// Build the matrix profile
void FELMRigidJoint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> lm;
	UnpackLM(lm);

	// add it to the pile
	M.build_add(lm);
}

void FELMRigidJoint::UnpackLM(vector<int>& lm)
{
	// get the displacement dofs
	FEModel& fem = *GetFEModel();
	DOFS& dofs = GetFEModel()->GetDOFS();
	int dofX = dofs.GetDOF("x");
	int dofY = dofs.GetDOF("y");
	int dofZ = dofs.GetDOF("z");

	// we need to couple the dofs of node A, B, and the LMs
	FEMesh& mesh = fem.GetMesh();

	// add the dofs of rigid body A
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
	lm.push_back(m_LM[0]);
	lm.push_back(m_LM[1]);
	lm.push_back(m_LM[2]);
}

void FELMRigidJoint::Update(const std::vector<double>& ui)
{
	m_Lm.x += ui[m_LM[0]];
	m_Lm.y += ui[m_LM[1]];
	m_Lm.z += ui[m_LM[2]];
}
