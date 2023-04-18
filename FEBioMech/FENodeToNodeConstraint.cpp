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
#include "FENodeToNodeConstraint.h"
#include <FECore/FEMesh.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FEGlobalMatrix.h>

BEGIN_FECORE_CLASS(FENodeToNodeConstraint, FENLConstraint)
	ADD_PARAMETER(m_a, FE_RANGE_GREATER(0), "a");
	ADD_PARAMETER(m_b, FE_RANGE_GREATER(0), "b");
END_FECORE_CLASS();

FENodeToNodeConstraint::FENodeToNodeConstraint(FEModel* fem) : FENLConstraint(fem)
{
	m_a = m_b = -1;
	m_Lm = vec3d(0, 0, 0);
}

// allocate equations
int FENodeToNodeConstraint::InitEquations(int neq)
{
	// we allocate three equations
	m_LM.resize(3);
	m_LM[0] = neq;
	m_LM[1] = neq+1;
	m_LM[2] = neq+2;
	return 3;
}

void FENodeToNodeConstraint::UnpackLM(vector<int>& lm)
{
	// get the displacement dofs
	int dofX = GetDOFIndex("x");
	int dofY = GetDOFIndex("y");
	int dofZ = GetDOFIndex("z");

	// we need to couple the dofs of node A, B, and the LMs
	FEMesh& mesh = GetMesh();

	// add the dofs of node A
	FENode& node_a = mesh.Node(m_a - 1);
	lm.push_back(node_a.m_ID[dofX]);
	lm.push_back(node_a.m_ID[dofY]);
	lm.push_back(node_a.m_ID[dofZ]);

	// add the dofs of node B
	FENode& node_b = mesh.Node(m_b - 1);
	lm.push_back(node_b.m_ID[dofX]);
	lm.push_back(node_b.m_ID[dofY]);
	lm.push_back(node_b.m_ID[dofZ]);

	// add the LM equations
	lm.push_back(m_LM[0]);
	lm.push_back(m_LM[1]);
	lm.push_back(m_LM[2]);
}

void FENodeToNodeConstraint::Update(const std::vector<double>& ui)
{
	m_Lm.x += ui[m_LM[0]];
	m_Lm.y += ui[m_LM[1]];
	m_Lm.z += ui[m_LM[2]];
}

// The LoadVector function evaluates the "forces" that contribute to the residual of the system
void FENodeToNodeConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEMesh& mesh = GetMesh();
	vec3d ra = mesh.Node(m_a - 1).m_rt;
	vec3d rb = mesh.Node(m_b - 1).m_rt;
	vec3d c = ra - rb;

	vector<double> fe(9, 0.0);
	fe[0] = -m_Lm.x;
	fe[1] = -m_Lm.y;
	fe[2] = -m_Lm.z;
	fe[3] =  m_Lm.x;
	fe[4] =  m_Lm.y;
	fe[5] =  m_Lm.z;
	fe[6] = -c.x;
	fe[7] = -c.y;
	fe[8] = -c.z;

	vector<int> lm;
	UnpackLM(lm);

	R.Assemble(lm, fe);
}

// Evaluates the contriubtion to the stiffness matrix
void FENodeToNodeConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
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
void FENodeToNodeConstraint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> lm;
	UnpackLM(lm);

	// add it to the pile
	M.build_add(lm);
}
