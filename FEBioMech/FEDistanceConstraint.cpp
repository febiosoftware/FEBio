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
#include "FEDistanceConstraint.h"
#include <FECore/FELinearSystem.h>
#include "FEBioMech.h"
#include <FECore/FEMesh.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEDistanceConstraint, FENLConstraint);
	ADD_PARAMETER(m_blaugon, "laugon" ); 
	ADD_PARAMETER(m_atol   , "augtol" );
	ADD_PARAMETER(m_eps    , "penalty");
	ADD_PARAMETER(m_node   , 2, "node");
	ADD_PARAMETER(m_nminaug, "minaug");
	ADD_PARAMETER(m_nmaxaug, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEDistanceConstraint::FEDistanceConstraint(FEModel* pfem) : FENLConstraint(pfem), m_dofU(pfem)
{
	m_eps = 0.0;
	m_atol = 0.01;
	m_blaugon = false;
	m_node[0] = -1;
	m_node[1] = -1;
	m_l0 = 0.0;
	m_Lm = 0.0;
	m_nminaug = 0;
	m_nmaxaug = 10;

	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
}

//-----------------------------------------------------------------------------
//! Initializes data structures. 
bool FEDistanceConstraint::Init()
{
	// get the FE mesh
	FEMesh& mesh = GetMesh();
	int NN = mesh.Nodes();

	// make sure the nodes are valid
	// (remember, they are one-based since they are defined in the input file)
	if ((m_node[0] <= 0)||(m_node[0] > NN)) return false;
	if ((m_node[1] <= 0)||(m_node[1] > NN)) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::Activate()
{
	// don't forget to call base class
	FENLConstraint::Activate();

	// get the FE mesh
	FEMesh& mesh = GetMesh();
	int NN = mesh.Nodes();

	// get the initial position of the two nodes
	vec3d ra = mesh.Node(m_node[0] - 1).m_rt;
	vec3d rb = mesh.Node(m_node[1] - 1).m_rt;

	// set the initial length
	m_l0 = (ra - rb).norm();
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	// get the FE mesh
	FEMesh& mesh = GetMesh();

	// get the two nodes
	FENode& nodea = mesh.Node(m_node[0] - 1);
	FENode& nodeb = mesh.Node(m_node[1] - 1);

	// get the current position of the two nodes
	vec3d ra = nodea.m_rt;
	vec3d rb = nodeb.m_rt;

	// calculate the force
	double lt = (ra - rb).norm();
	double Lm = m_Lm + m_eps*(lt - m_l0);
	vec3d Fc = (ra - rb)*(Lm/lt);

	// setup the "element" force vector
	vector<double> fe(6);
	fe[0] = -Fc.x;
	fe[1] = -Fc.y;
	fe[2] = -Fc.z;
	fe[3] =  Fc.x;
	fe[4] =  Fc.y;
	fe[5] =  Fc.z;

	// setup the LM vector
	vector<int> lm(6);
	lm[0] = nodea.m_ID[m_dofU[0]];
	lm[1] = nodea.m_ID[m_dofU[1]];
	lm[2] = nodea.m_ID[m_dofU[2]];
	lm[3] = nodeb.m_ID[m_dofU[0]];
	lm[4] = nodeb.m_ID[m_dofU[1]];
	lm[5] = nodeb.m_ID[m_dofU[2]];

	// setup element vector
	vector<int> en(2);
	en[0] = m_node[0] - 1;
	en[1] = m_node[1] - 1;

	// add element force vector to global force vector
	R.Assemble(en, lm, fe);
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	// get the FE mesh
	FEMesh& mesh = GetMesh();

	// get the two nodes
	FENode& nodea = mesh.Node(m_node[0] - 1);
	FENode& nodeb = mesh.Node(m_node[1] - 1);

	// get the current position of the two nodes
	vec3d ra = nodea.m_rt;
	vec3d rb = nodeb.m_rt;
	vec3d rab = ra - rb;

	// calculate the Lagrange mulitplier
	double lt = rab.norm();
	double l0 = m_l0;
	double l3 = lt*lt*lt;
	double Lm = m_Lm + m_eps*(lt - l0);

	// calculate the stiffness
	mat3d kab;
	kab[0][0] = (m_eps*l0*rab.x*rab.x/l3 + Lm/lt);
	kab[1][1] = (m_eps*l0*rab.y*rab.y/l3 + Lm/lt);
	kab[2][2] = (m_eps*l0*rab.z*rab.z/l3 + Lm/lt);
	kab[0][1] = kab[1][0] = (m_eps*l0*rab.x*rab.y/l3);
	kab[0][2] = kab[2][0] = (m_eps*l0*rab.x*rab.z/l3);
	kab[1][2] = kab[2][1] = (m_eps*l0*rab.y*rab.z/l3);

	// element stiffness matrix
	FEElementMatrix ke;
	ke.resize(6, 6);
	ke.zero();
	ke[0][0] = kab[0][0]; ke[0][1] = kab[0][1]; ke[0][2] = kab[0][2]; ke[0][3] = -kab[0][0]; ke[0][4] = -kab[0][1]; ke[0][5] = -kab[0][2];
	ke[1][0] = kab[1][0]; ke[1][1] = kab[1][1]; ke[1][2] = kab[1][2]; ke[1][3] = -kab[1][0]; ke[1][4] = -kab[1][1]; ke[1][5] = -kab[1][2];
	ke[2][0] = kab[2][0]; ke[2][1] = kab[2][1]; ke[2][2] = kab[2][2]; ke[2][3] = -kab[2][0]; ke[2][4] = -kab[2][1]; ke[2][5] = -kab[2][2];

	ke[3][0] = -kab[0][0]; ke[3][1] = -kab[0][1]; ke[3][2] = -kab[0][2]; ke[3][3] = kab[0][0]; ke[3][4] = kab[0][1]; ke[3][5] = kab[0][2];
	ke[4][0] = -kab[1][0]; ke[4][1] = -kab[1][1]; ke[4][2] = -kab[1][2]; ke[4][3] = kab[1][0]; ke[4][4] = kab[1][1]; ke[4][5] = kab[1][2];
	ke[5][0] = -kab[2][0]; ke[5][1] = -kab[2][1]; ke[5][2] = -kab[2][2]; ke[5][3] = kab[2][0]; ke[5][4] = kab[2][1]; ke[5][5] = kab[2][2];

	// setup the LM vector
	vector<int> lm(6);
	lm[0] = nodea.m_ID[m_dofU[0]];
	lm[1] = nodea.m_ID[m_dofU[1]];
	lm[2] = nodea.m_ID[m_dofU[2]];
	lm[3] = nodeb.m_ID[m_dofU[0]];
	lm[4] = nodeb.m_ID[m_dofU[1]];
	lm[5] = nodeb.m_ID[m_dofU[2]];

	// setup element vector
	vector<int> en(2);
	en[0] = m_node[0] - 1;
	en[1] = m_node[1] - 1;

	// assemble element matrix in global stiffness matrix
	ke.SetNodes(en);
	ke.SetIndices(lm);
	LS.Assemble(ke);
}

//-----------------------------------------------------------------------------
bool FEDistanceConstraint::Augment(int naug, const FETimeInfo& tp)
{
	// make sure we are augmenting
	if ((m_blaugon == false) || (m_atol <= 0.0)) return true;

	// get the FE mesh
	FEMesh& mesh = GetMesh();

	// get the two nodes
	FENode& nodea = mesh.Node(m_node[0] - 1);
	FENode& nodeb = mesh.Node(m_node[1] - 1);

	// get the current position of the two nodes
	vec3d ra = nodea.m_rt;
	vec3d rb = nodeb.m_rt;

	// calculate the Lagrange multipler
	double l = (ra - rb).norm();
	double Lm = m_Lm + m_eps*(l - m_l0);

	// calculate force
	vec3d Fc = (ra - rb)*(Lm/l);

	// calculate relative error
	bool bconv = false;
	double err = fabs((Lm - m_Lm)/Lm);
	if (err < m_atol) bconv = true;

	feLog("\ndistance constraint:\n");
	feLog("\tmultiplier= %lg (error = %lg / %lg)\n", Lm, err, m_atol);
	feLog("\tforce     = %lg, %lg, %lg\n", Fc.x, Fc.y, Fc.z);
	feLog("\tdistance  = %lg (L0 = %lg)\n", l, m_l0);

	// check convergence
	if (m_nminaug > naug) bconv = false;
	if (m_nmaxaug <= naug) bconv = true;

	// update Lagrange multiplier
	// (only when we did not converge)
	if (bconv == false) m_Lm = Lm;

	return bconv;
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	FEMesh& mesh = GetMesh();
	vector<int> lm(6);
	FENode& n0 = mesh.Node(m_node[0] - 1);
	lm[0] = n0.m_ID[m_dofU[0]];
	lm[1] = n0.m_ID[m_dofU[1]];
	lm[2] = n0.m_ID[m_dofU[2]];
	FENode& n1 = mesh.Node(m_node[1] - 1);
	lm[3] = n1.m_ID[m_dofU[0]];
	lm[4] = n1.m_ID[m_dofU[1]];
	lm[5] = n1.m_ID[m_dofU[2]];
    M.build_add(lm);
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::Serialize(DumpStream& ar)
{
	FENLConstraint::Serialize(ar);
	ar & m_Lm; 
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::Reset()
{
	m_Lm = 0.0;
}
