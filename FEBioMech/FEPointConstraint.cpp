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
#include "FEPointConstraint.h"
#include "FECore/FEMesh.h"
#include "FEBioMech.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPointConstraint, FENLConstraint)
	ADD_PARAMETER(m_eps    , "penalty");
	ADD_PARAMETER(m_node_id, "node"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPointConstraint::FEPointConstraint(FEModel* pfem) : FENLConstraint(pfem), m_dofU(pfem)
{
	m_node_id = -1;
	m_eps = 0.0;

	m_node = -1;
	m_pel = 0;

	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
}

//-----------------------------------------------------------------------------
bool FEPointConstraint::Init()
{
	FEMesh& m = GetMesh();
	if ((m_node_id <= 0)||(m_node_id > m.Nodes())) return false;

	// get the nodal position in the reference state
	m_node = m_node_id - 1;
	vec3d r = m.Node(m_node).m_r0;

	// find the element in which this node lies
	m_pel = m.FindSolidElement(r, m_rs);
	if (m_pel == 0) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEPointConstraint::Serialize(DumpStream& ar)
{
	FENLConstraint::Serialize(ar);
	// TODO: implement this
}

//-----------------------------------------------------------------------------
void FEPointConstraint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	FEMesh& mesh = GetMesh();
	vector<int> lm(3*9);
	FENode& n0 = mesh.Node(m_node);
	lm[0] = n0.m_ID[m_dofU[0]];
	lm[1] = n0.m_ID[m_dofU[1]];
	lm[2] = n0.m_ID[m_dofU[2]];
	for (int i=0; i<8; ++i)
	{
		FENode& ni = mesh.Node(m_pel->m_node[i]);
		lm[3*(i+1)  ] = ni.m_ID[m_dofU[0]];
		lm[3*(i+1)+1] = ni.m_ID[m_dofU[1]];
		lm[3*(i+1)+2] = ni.m_ID[m_dofU[2]];
	}
	M.build_add(lm);
}

//-----------------------------------------------------------------------------
void FEPointConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	int i;
	FEMesh& m = GetMesh();

	// calculate H matrix
	double H[9], *r = m_rs;
	H[0] = 1.0;
	H[1] = -0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
	H[2] = -0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
	H[3] = -0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
	H[4] = -0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
	H[5] = -0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
	H[6] = -0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
	H[7] = -0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
	H[8] = -0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);

	// get the nodal position
	vec3d x[9];
	x[0] = m.Node(m_node).m_rt;
	for (i=0; i<8; ++i) x[i+1] = m.Node(m_pel->m_node[i]).m_rt;

	// calculate the constraint
	vec3d c(0,0,0);
	for (i=0; i<9; ++i) c += x[i]*H[i];

	// calculate the force
	vec3d T = c*m_eps;

	// setup the LM matrix
	vector<int> LM(3*9), en(9);
	en[0] = m_node;
	LM[0] = m.Node(m_node).m_ID[m_dofU[0]];
	LM[1] = m.Node(m_node).m_ID[m_dofU[1]];
	LM[2] = m.Node(m_node).m_ID[m_dofU[2]];
	for (i=0; i<8; ++i)
	{
		en[i+1] = m_pel->m_node[i];
		FENode& node = m.Node(en[i+1]);
		LM[(i+1)*3  ] = node.m_ID[m_dofU[0]];
		LM[(i+1)*3+1] = node.m_ID[m_dofU[1]];
		LM[(i+1)*3+2] = node.m_ID[m_dofU[2]];
	}

	// set up nodal force vector
	vector<double> fe(3*9);
	for (int i=0; i<9; ++i)
	{
		fe[3*i  ] = -T.x*H[i];
		fe[3*i+1] = -T.y*H[i];
		fe[3*i+2] = -T.z*H[i];
	}

	// assemble residual
	R.Assemble(en, LM, fe);
}

//-----------------------------------------------------------------------------
void FEPointConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEMesh& m = GetMesh();

	// calculate H matrix
	double H[9], *r = m_rs;
	H[0] = 1.0;
	H[1] = -0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
	H[2] = -0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
	H[3] = -0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
	H[4] = -0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
	H[5] = -0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
	H[6] = -0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
	H[7] = -0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
	H[8] = -0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);


	// setup the LM matrix
	vector<int> LM(3*9), en(9);
	en[0] = m_node;
	LM[0] = m.Node(m_node).m_ID[m_dofU[0]];
	LM[1] = m.Node(m_node).m_ID[m_dofU[1]];
	LM[2] = m.Node(m_node).m_ID[m_dofU[2]];
	for (int i=0; i<8; ++i)
	{
		en[i+1] = m_pel->m_node[i];
		FENode& node = m.Node(en[i+1]);
		LM[(i+1)*3  ] = node.m_ID[m_dofU[0]];
		LM[(i+1)*3+1] = node.m_ID[m_dofU[1]];
		LM[(i+1)*3+2] = node.m_ID[m_dofU[2]];
	}

	// setup stiffness matrix
	int ndof = 3*9;
	FEElementMatrix ke(ndof, ndof); ke.zero();
	for (int i=0; i<9; ++i)
		for (int j=0; j<9; ++j)
		{
			ke[3*i  ][3*j  ] = m_eps*H[i]*H[j];
			ke[3*i+1][3*j+1] = m_eps*H[i]*H[j];
			ke[3*i+2][3*j+2] = m_eps*H[i]*H[j];
		}
	ke.SetIndices(LM);
	ke.SetNodes(en);

	// assemble stiffness matrix
	LS.Assemble(ke);
}
