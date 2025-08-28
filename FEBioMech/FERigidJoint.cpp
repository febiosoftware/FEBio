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
	ADD_PARAMETER(m_laugon , "laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0LAGMULT\0");
	ADD_PARAMETER(m_atol   , "tolerance");
	ADD_PARAMETER(m_eps    , "penalty"  );
	ADD_PARAMETER(m_q0     , "joint"    );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidJoint::FERigidJoint(FEModel* pfem) : FERigidConnector(pfem)
{
	static int count = 1;
	m_nID = count++;

	m_laugon = FECore::AUGLAG_METHOD;		// Augmented Lagrangian by default for backward compatibility
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

	// we make sure we have a non-zero penalty for penalty and auglag method
	if ((m_laugon != FECore::LAGMULT_METHOD) && (m_eps == 0.0)) return false;

	return true;
}

//-----------------------------------------------------------------------------
// allocate equations
int FERigidJoint::InitEquations(int neq)
{
	m_LM.resize(3, -1);
	if (m_laugon == FECore::LAGMULT_METHOD)
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
	double alpha = tp.alpha;

	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	// body A
	vec3d rat = RBa.m_rt;
	vec3d rap = RBa.m_rp;
	vec3d ra = rat * alpha + rap*(1 - alpha);

	vec3d zat = RBa.GetRotation() * m_qa0;
	vec3d zap = RBa.m_qp * m_qa0;
	vec3d za = zat * alpha + zap * (1 - alpha);

	// body b
	vec3d rbt = RBb.m_rt;
	vec3d rbp = RBb.m_rp;
	vec3d rb = rbt*alpha + rbp*(1 - alpha);

	vec3d zbt = RBb.GetRotation() * m_qb0;
	vec3d zbp = RBb.m_qp * m_qb0;
	vec3d zb = zbt * alpha + zbp * (1 - alpha);

	// constraint
	vec3d c = ra + za - rb - zb;

	// forces
	vec3d F = m_L + c * m_eps;
	if (m_laugon == FECore::LAGMULT_METHOD) F = m_F;

	vec3d Fa =  F;
	vec3d Fb = -F;
	vec3d Ma =  (za ^ F);
	vec3d Mb = -(zb ^ F);

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

	double alpha = tp.alpha;

	// body A
	vec3d rat = RBa.m_rt;
	vec3d rap = RBa.m_rp;
	vec3d ra = rat * alpha + rap * (1 - alpha);

	vec3d zat = RBa.GetRotation() * m_qa0;
	vec3d zap = RBa.m_qp * m_qa0;
	vec3d za = zat * alpha + zap * (1 - alpha);

	// body b
	vec3d rbt = RBb.m_rt;
	vec3d rbp = RBb.m_rp;
	vec3d rb = rbt * alpha + rbp * (1 - alpha);

	vec3d zbt = RBb.GetRotation() * m_qb0;
	vec3d zbp = RBb.m_qp * m_qb0;
	vec3d zb = zbt * alpha + zbp * (1 - alpha);

	FEElementMatrix ke;

	if (m_laugon != FECore::LAGMULT_METHOD)
	{
		ke.resize(12, 12);
		ke.zero();

		// constraint
		vec3d c = ra + za - rb - zb;

		// forces
		vec3d F = m_L + c * m_eps;
		mat3da Fhat(F);

		mat3da zahat(za);
		mat3da zathat(zat);

		mat3da zbhat(zb);
		mat3da zbthat(zbt);

		mat3dd I(1.0);

		ke.set(0, 0,      I); ke.set(0, 3,       -zathat); ke.set(0, 6,     -I); ke.set(0, 9,        zbthat);
		ke.set(3, 0,  zahat); ke.set(3, 3, -zahat*zathat); ke.set(3, 6, -zahat); ke.set(3, 9,  zahat*zbthat);
		ke.set(6, 0,     -I); ke.set(6, 3,        zathat); ke.set(6, 6,      I); ke.set(6, 9,       -zbthat);
		ke.set(9, 0, -zbhat); ke.set(9, 3,  zbhat*zathat); ke.set(9, 6,  zbhat); ke.set(9, 9, -zbhat*zbthat);

		// scale by penalty factor
		ke *= m_eps;

		ke.add(3, 3,  Fhat*zathat);
		ke.add(9, 9, -Fhat*zbthat);

		ke *= alpha;
	}
	else
	{
		ke.resize(15, 15);
		ke.zero();

		mat3dd I(1.0);
		mat3da yaT(za);
		mat3da ybT(zb);

		mat3da Fhat(m_F);

		ke.add_symm(0, 12,   I);
		ke.add_symm(3, 12,  yaT);
		ke.add_symm(6, 12,  -I);
		ke.add_symm(9, 12, -ybT);

		ke.add(3, 3,  Fhat*yaT);
		ke.add(9, 9, -Fhat*ybT);
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
	if (m_laugon == FECore::PENALTY_METHOD) return true;

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
	if (m_laugon == FECore::LAGMULT_METHOD)
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
	if (m_laugon != FECore::LAGMULT_METHOD)
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
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		lm.push_back(m_LM[0]);
		lm.push_back(m_LM[1]);
		lm.push_back(m_LM[2]);
	}
}

//-----------------------------------------------------------------------------
void FERigidJoint::Update(const std::vector<double>& ui)
{
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		m_F.x += ui[m_LM[0]];
		m_F.y += ui[m_LM[1]];
		m_F.z += ui[m_LM[2]];
	}
}
