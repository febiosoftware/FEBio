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
END_FECORE_CLASS();

FELMRigidJoint::FELMRigidJoint(FEModel* fem) : FELMConstraint(fem)
{
	m_nRBa = -1;
	m_nRBb = -1;
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
	vec3d c = ra - rb;

	vector<double> fe(9, 0.0);
	fe[0] = -m_Lm.x;
	fe[1] = -m_Lm.y;
	fe[2] = -m_Lm.z;
	fe[3] = m_Lm.x;
	fe[4] = m_Lm.y;
	fe[5] = m_Lm.z;
	fe[6] = -c.x;
	fe[7] = -c.y;
	fe[8] = -c.z;

	vector<int> lm;
	UnpackLM(lm);

	R.Assemble(lm, fe);
}

// Evaluates the contriubtion to the stiffness matrix
void FELMRigidJoint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEElementMatrix ke(9, 9);
	ke.zero();

	ke[0][6] = ke[6][0] = 1;
	ke[1][7] = ke[7][1] = 1;
	ke[2][8] = ke[8][2] = 1;

	ke[3][6] = ke[6][3] = -1;
	ke[4][7] = ke[7][4] = -1;
	ke[5][8] = ke[8][5] = -1;

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

	// add the dofs of node A
	lm.push_back(m_rbA->m_LM[0]);
	lm.push_back(m_rbA->m_LM[1]);
	lm.push_back(m_rbA->m_LM[2]);

	// add the dofs of node B
	lm.push_back(m_rbB->m_LM[0]);
	lm.push_back(m_rbB->m_LM[1]);
	lm.push_back(m_rbB->m_LM[2]);

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
