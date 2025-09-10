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
#include "FERigidLock.h"
#include "FERigidBody.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidLock, FERigidConnector);
	ADD_PARAMETER(m_laugon  , "laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0LAGMULT\0");
	ADD_PARAMETER(m_atol    , "tolerance"     );
	ADD_PARAMETER(m_gtol    , "gaptol"        );
	ADD_PARAMETER(m_qtol    , "angtol"        );
	ADD_PARAMETER(m_eps     , "force_penalty" );
	ADD_PARAMETER(m_ups     , "moment_penalty");
	ADD_PARAMETER(m_q0      , "joint_origin"  );
	ADD_PARAMETER(m_e0[0]   , "first_axis"    );
	ADD_PARAMETER(m_e0[1]   , "second_axis"   );
	ADD_PARAMETER(m_naugmin , "minaug"        );
	ADD_PARAMETER(m_naugmax , "maxaug"        );
	ADD_PARAMETER(m_bautopen, "auto_penalty");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidLock::FERigidLock(FEModel* pfem) : FERigidConnector(pfem)
{
    m_eps = 1.0;
    m_ups = 1.0;

    m_nID = m_ncount++;
	m_laugon = FECore::AUGLAG_METHOD; // for backward compatibility
    m_atol = 0;
    m_gtol = 0;
    m_qtol = 0;
    m_naugmin = 0;
    m_naugmax = 10;
	m_bautopen = false;
    m_eps = m_ups = 1.0;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidLock::Init()
{
    // initialize joint basis
    m_e0[0].unit();
    m_e0[2] = m_e0[0] ^ m_e0[1]; m_e0[2].unit();
    m_e0[1] = m_e0[2] ^ m_e0[0]; m_e0[1].unit();
    
    // reset force
    m_F = vec3d(0,0,0); m_L = vec3d(0,0,0);
    m_M = vec3d(0,0,0); m_U = vec3d(0,0,0);
	m_U1 = m_U2 = vec3d(0, 0, 0);
    
	// base class first
	if (FERigidConnector::Init() == false) return false;

    m_qa0 = m_q0 - m_rbA->m_r0;
    m_qb0 = m_q0 - m_rbB->m_r0;
    
    m_ea0[0] = m_e0[0]; m_ea0[1] = m_e0[1]; m_ea0[2] = m_e0[2];
    m_eb0[0] = m_e0[0]; m_eb0[1] = m_e0[1]; m_eb0[2] = m_e0[2];
    
    return true;
}

// allocate equations
int FERigidLock::InitEquations(int neq)
{
	const int LMeq = 9;
	m_LM.resize(LMeq, -1);
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		for (int i=0; i< LMeq; ++i) m_LM[i] = neq + i;
		return LMeq;
	}
	else return 0;
}

// Build the matrix profile
void FERigidLock::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> lm;
	UnpackLM(lm);

	// add it to the pile
	M.build_add(lm);
}

void FERigidLock::UnpackLM(vector<int>& lm)
{
	// add the dofs of rigid body A
	lm.resize(21);
	lm[0] = m_rbA->m_LM[0];
	lm[1] = m_rbA->m_LM[1];
	lm[2] = m_rbA->m_LM[2];
	lm[3] = m_rbA->m_LM[3];
	lm[4] = m_rbA->m_LM[4];
	lm[5] = m_rbA->m_LM[5];

	// add the dofs of rigid body B
	lm[ 6] = m_rbB->m_LM[0];
	lm[ 7] = m_rbB->m_LM[1];
	lm[ 8] = m_rbB->m_LM[2];
	lm[ 9] = m_rbB->m_LM[3];
	lm[10] = m_rbB->m_LM[4];
	lm[11] = m_rbB->m_LM[5];

	// add the LM equations
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		for (int i = 0; i < m_LM.size(); ++i) lm[12 + i] = m_LM[i];
	}
}

//-----------------------------------------------------------------------------
void FERigidLock::Serialize(DumpStream& ar)
{
    FERigidConnector::Serialize(ar);
    ar & m_qa0 & m_qb0;
    ar & m_L & m_U & m_U1 & m_U2;
    ar & m_e0;
    ar & m_ea0;
    ar & m_eb0;
	ar & m_LM;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidLock::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<double> fa(6);
    vector<double> fb(6);
    
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    double alpha = tp.alphaf;
    
    // body A
    vec3d ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
	eat[0] = m_ea0[0]; RBa.GetRotation().RotateVector(eat[0]);
	eat[1] = m_ea0[1]; RBa.GetRotation().RotateVector(eat[1]);
	eat[2] = m_ea0[2]; RBa.GetRotation().RotateVector(eat[2]);
    eap[0] = m_ea0[0]; RBa.m_qp.RotateVector(eap[0]);
    eap[1] = m_ea0[1]; RBa.m_qp.RotateVector(eap[1]);
    eap[2] = m_ea0[2]; RBa.m_qp.RotateVector(eap[2]);
    ea[0] = eat[0]*alpha + eap[0]*(1-alpha);
    ea[1] = eat[1]*alpha + eap[1]*(1-alpha);
    ea[2] = eat[2]*alpha + eap[2]*(1-alpha);
    
    // body b
    vec3d rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
	ebt[0] = m_eb0[0]; RBb.GetRotation().RotateVector(ebt[0]);
	ebt[1] = m_eb0[1]; RBb.GetRotation().RotateVector(ebt[1]);
	ebt[2] = m_eb0[2]; RBb.GetRotation().RotateVector(ebt[2]);
    ebp[0] = m_eb0[0]; RBb.m_qp.RotateVector(ebp[0]);
    ebp[1] = m_eb0[1]; RBb.m_qp.RotateVector(ebp[1]);
    ebp[2] = m_eb0[2]; RBb.m_qp.RotateVector(ebp[2]);
    eb[0] = ebt[0]*alpha + ebp[0]*(1-alpha);
    eb[1] = ebt[1]*alpha + ebp[1]*(1-alpha);
    eb[2] = ebt[2]*alpha + ebp[2]*(1-alpha);
    
    // incremental compound rotation of B w.r.t. A
	if (m_laugon != FECore::LAGMULT_METHOD)
	{
		vec3d vth = ((ea[0] ^ eb[0]) + (ea[1] ^ eb[1]) + (ea[2] ^ eb[2])) / 2;

		vec3d c = rb + zb - ra - za;
		if (m_laugon != FECore::LAGMULT_METHOD) m_F = m_L + c * m_eps;

		vec3d ksi = vth;
		if (m_laugon != FECore::LAGMULT_METHOD) m_M = m_U + ksi * m_ups;

		fa[0] = m_F.x;
		fa[1] = m_F.y;
		fa[2] = m_F.z;

		fa[3] = za.y * m_F.z - za.z * m_F.y + m_M.x;
		fa[4] = za.z * m_F.x - za.x * m_F.z + m_M.y;
		fa[5] = za.x * m_F.y - za.y * m_F.x + m_M.z;

		fb[0] = -m_F.x;
		fb[1] = -m_F.y;
		fb[2] = -m_F.z;

		fb[3] = -zb.y * m_F.z + zb.z * m_F.y - m_M.x;
		fb[4] = -zb.z * m_F.x + zb.x * m_F.z - m_M.y;
		fb[5] = -zb.x * m_F.y + zb.y * m_F.x - m_M.z;

		for (int i = 0; i < 6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
		for (int i = 0; i < 6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];
	}
	else
	{
		vec3d Ma = (za ^ m_F) + (ea[0] ^ m_U1) + (ea[1] ^ m_U2);
		vec3d Mb = (zb ^ m_F) + (eb[0] ^ m_U1) + (eb[1] ^ m_U2);

		fa[0] = -m_F.x;
		fa[1] = -m_F.y;
		fa[2] = -m_F.z;
		fa[3] = -Ma.x;
		fa[4] = -Ma.y;
		fa[5] = -Ma.z;

		fb[0] =  m_F.x;
		fb[1] =  m_F.y;
		fb[2] =  m_F.z;
		fb[3] =  Mb.x;
		fb[4] =  Mb.y;
		fb[5] =  Mb.z;

		for (int i = 0; i < 6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
		for (int i = 0; i < 6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];

		// translational constraint
		vec3d c = rb + zb - ra - za;
		R[m_LM[0]] += c.x;
		R[m_LM[1]] += c.y;
		R[m_LM[2]] += c.z;

		// rotational constraint
		vec3d ksi1 = eb[0] - ea[0];
		vec3d ksi2 = eb[1] - ea[1];
		R[m_LM[3]] += ksi1.x;
		R[m_LM[4]] += ksi1.y;
		R[m_LM[5]] += ksi1.z;
		R[m_LM[6]] += ksi2.x;
		R[m_LM[7]] += ksi2.y;
		R[m_LM[8]] += ksi2.z;
	}
    
    RBa.m_Fr -= vec3d(fa[0],fa[1],fa[2]);
    RBa.m_Mr -= vec3d(fa[3],fa[4],fa[5]);
    RBb.m_Fr -= vec3d(fb[0],fb[1],fb[2]);
    RBb.m_Mr -= vec3d(fb[3],fb[4],fb[5]);
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidLock::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    double alpha = tp.alphaf;
    
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    // body A
	quatd Qat = RBa.GetRotation();
	quatd Qap = RBa.m_qp;
	vec3d rat = RBa.m_rt;
	vec3d rap = RBa.m_rp;
    vec3d ra = rat*alpha + rap*(1-alpha);
	vec3d zat = Qat*m_qa0;
    vec3d zap = Qap*m_qa0;
    vec3d za = zat*alpha + zap*(1-alpha);
	eat[0] = m_ea0[0]; RBa.GetRotation().RotateVector(eat[0]);
	eat[1] = m_ea0[1]; RBa.GetRotation().RotateVector(eat[1]);
	eat[2] = m_ea0[2]; RBa.GetRotation().RotateVector(eat[2]);
    eap[0] = m_ea0[0]; RBa.m_qp.RotateVector(eap[0]);
    eap[1] = m_ea0[1]; RBa.m_qp.RotateVector(eap[1]);
    eap[2] = m_ea0[2]; RBa.m_qp.RotateVector(eap[2]);
    ea[0] = eat[0]*alpha + eap[0]*(1-alpha);
    ea[1] = eat[1]*alpha + eap[1]*(1-alpha);
    ea[2] = eat[2]*alpha + eap[2]*(1-alpha);
    mat3d zahat; zahat.skew(za);
    mat3d zathat; zathat.skew(zat);
    
    // body b
	quatd Qbt = RBb.GetRotation();
	quatd Qbp = RBb.m_qp;
	vec3d rbt = RBb.m_rt;
	vec3d rbp = RBb.m_rp;
	vec3d rb = rbt*alpha + rbp*(1-alpha);
	vec3d zbt = Qbt*m_qb0;
    vec3d zbp = Qbp*m_qb0;
    vec3d zb = zbt*alpha + zbp*(1-alpha);
	ebt[0] = m_eb0[0]; RBb.GetRotation().RotateVector(ebt[0]);
	ebt[1] = m_eb0[1]; RBb.GetRotation().RotateVector(ebt[1]);
	ebt[2] = m_eb0[2]; RBb.GetRotation().RotateVector(ebt[2]);
    ebp[0] = m_eb0[0]; RBb.m_qp.RotateVector(ebp[0]);
    ebp[1] = m_eb0[1]; RBb.m_qp.RotateVector(ebp[1]);
    ebp[2] = m_eb0[2]; RBb.m_qp.RotateVector(ebp[2]);
    eb[0] = ebt[0]*alpha + ebp[0]*(1-alpha);
    eb[1] = ebt[1]*alpha + ebp[1]*(1-alpha);
    eb[2] = ebt[2]*alpha + ebp[2]*(1-alpha);
    mat3d zbhat; zbhat.skew(zb);
    mat3d zbthat; zbthat.skew(zbt);
    
    mat3d eahat[3], ebhat[3], eathat[3], ebthat[3];
    for (int j=0; j<3; ++j) {
        eahat[j] = skew(ea[j]);
        ebhat[j] = skew(eb[j]);
        eathat[j] = skew(eat[j]);
        ebthat[j] = skew(ebt[j]);
    }

	if (m_laugon != FECore::LAGMULT_METHOD)
	{
		vector<int> LM(12);
		FEElementMatrix ke; ke.resize(12, 12);
		ke.zero();

		vec3d c = rb + zb - ra - za;
		m_F = m_L + c * m_eps;

		// incremental compound rotation of B w.r.t. A
		vec3d vth = ((ea[0] ^ eb[0]) + (ea[1] ^ eb[1]) + (ea[2] ^ eb[2])) / 2;
		vec3d ksi = vth;
		m_M = m_U + ksi * m_ups;

		mat3da Fhat(m_F);

		mat3dd I(1);
		mat3d Wba = (ebhat[0] * eathat[0] + ebhat[1] * eathat[1] + ebhat[2] * eathat[2]) / 2;
		mat3d Wab = (eahat[0] * ebthat[0] + eahat[1] * ebthat[1] + eahat[2] * ebthat[2]) / 2;

		// row 1
		ke.set(0, 0, I * (-m_eps));
		ke.set(0, 3, zathat * (m_eps));
		ke.set(0, 6, I * (m_eps));
		ke.set(0, 9, zbthat * (-m_eps));

		// row 2
		ke.set(3, 0, zahat * (-m_eps));
		ke.set(3, 3, zahat * zathat * (m_eps)+Fhat * zathat + Wba * m_ups);
		ke.set(3, 6, zahat * (m_eps));
		ke.set(3, 9, zahat * zbthat * (-m_eps) + Wab * (-m_ups));

		// row 3
		ke.set(6, 0, I * (m_eps));
		ke.set(6, 3, zathat * (-m_eps));
		ke.set(6, 6, I * (-m_eps));
		ke.set(6, 9, zbthat * (m_eps));

		// row 4
		ke.set(9, 0, zbhat * (m_eps));
		ke.set(9, 3, zbhat * zathat * (-m_eps) + Wba * (-m_ups));
		ke.set(9, 6, zbhat * (-m_eps));
		ke.set(9, 9, zbhat * zbthat * (m_eps)-Fhat * zbthat + Wab * m_ups);

		ke *= -alpha;

		for (int j = 0; j < 6; ++j)
		{
			LM[j] = RBa.m_LM[j];
			LM[j + 6] = RBb.m_LM[j];
		}

		ke.SetIndices(LM);
		LS.Assemble(ke);
	}
	else
	{
		FEElementMatrix ke; ke.resize(21, 21);
		ke.zero();

		mat3dd I(1.0);

		// translation constraint
		mat3da Fhat(m_F);
		ke.sub(3, 3, -Fhat * zahat);
		ke.sub(9, 9,  Fhat * zbhat);
		ke.sub(0, 12, -I);
		ke.sub(3, 12, -zahat);
		ke.sub(6, 12,  I);
		ke.sub(9, 12,  zbhat);
		ke.sub(12, 0, -I);
		ke.sub(12, 3, zahat);
		ke.sub(12, 6,  I);
		ke.sub(12, 9, -zbhat);

		// rotational constraint 1
		mat3da U1hat(m_U1);
		ke.sub(3, 3, -U1hat*eahat[0]);
		ke.sub(9, 9,  U1hat*ebhat[0]);
		ke.sub(3, 15, -eahat[0]);
		ke.sub(9, 15,  ebhat[0]);
		ke.sub(15, 3,  eahat[0]);
		ke.sub(15, 9, -ebhat[0]);

		// rotational constraint 2
		mat3da U2hat(m_U2);
		ke.sub(3, 3, -U2hat*eahat[1]);
		ke.sub(9, 9,  U2hat*ebhat[1]);
		ke.sub(3, 18, -eahat[1]);
		ke.sub(9, 18,  ebhat[1]);
		ke.sub(18, 3,  eahat[1]);
		ke.sub(18, 9, -ebhat[1]);

		vector<int> LM;
		UnpackLM(LM);
		ke.SetIndices(LM);
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
bool FERigidLock::Augment(int naug, const FETimeInfo& tp)
{
	if (m_laugon != FECore::AUGLAG_METHOD) return true;

    vec3d ra, rb, qa, qb, c, ksi, Lm;
    vec3d za, zb;
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    double normF0, normF1;
    vec3d Um;
    double normM0, normM1;
    bool bconv = true;
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    double alpha = tp.alphaf;
    
    ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*alpha + zap*(1-alpha);
	eat[0] = m_ea0[0]; RBa.GetRotation().RotateVector(eat[0]);
	eat[1] = m_ea0[1]; RBa.GetRotation().RotateVector(eat[1]);
	eat[2] = m_ea0[2]; RBa.GetRotation().RotateVector(eat[2]);
    eap[0] = m_ea0[0]; RBa.m_qp.RotateVector(eap[0]);
    eap[1] = m_ea0[1]; RBa.m_qp.RotateVector(eap[1]);
    eap[2] = m_ea0[2]; RBa.m_qp.RotateVector(eap[2]);
    ea[0] = eat[0]*alpha + eap[0]*(1-alpha);
    ea[1] = eat[1]*alpha + eap[1]*(1-alpha);
    ea[2] = eat[2]*alpha + eap[2]*(1-alpha);
    
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*alpha + zbp*(1-alpha);
	ebt[0] = m_eb0[0]; RBb.GetRotation().RotateVector(ebt[0]);
	ebt[1] = m_eb0[1]; RBb.GetRotation().RotateVector(ebt[1]);
	ebt[2] = m_eb0[2]; RBb.GetRotation().RotateVector(ebt[2]);
    ebp[0] = m_eb0[0]; RBb.m_qp.RotateVector(ebp[0]);
    ebp[1] = m_eb0[1]; RBb.m_qp.RotateVector(ebp[1]);
    ebp[2] = m_eb0[2]; RBb.m_qp.RotateVector(ebp[2]);
    eb[0] = ebt[0]*alpha + ebp[0]*(1-alpha);
    eb[1] = ebt[1]*alpha + ebp[1]*(1-alpha);
    eb[2] = ebt[2]*alpha + ebp[2]*(1-alpha);
    
    // incremental compound rotation of B w.r.t. A
    vec3d vth = ((ea[0] ^ eb[0]) + (ea[1] ^ eb[1]) + (ea[2] ^ eb[2]))/2;
    
    c = rb + zb - ra - za;
    
    normF0 = sqrt(m_L*m_L);
    
    // calculate trial multiplier
    Lm = m_L + c*m_eps;
    
    normF1 = sqrt(Lm*Lm);
    
    ksi = vth;
    
    normM0 = sqrt(m_U*m_U);
    
    // calculate trial multiplier
    Um = m_U + ksi*m_ups;
    
    normM1 = sqrt(Um*Um);
    
    // check convergence of constraints
    feLog(" rigid connector # %d (lock)\n", m_nID+1);
    feLog("                  CURRENT        REQUIRED\n");
    double pctn = 0;
    double gap = c.norm();
    double qap = ksi.norm();
    if (fabs(normF1) > 1e-10) pctn = fabs((normF1 - normF0)/normF1);
    if (m_atol) feLog("    force : %15le %15le\n", pctn, m_atol);
    else        feLog("    force : %15le        ***\n", pctn);
    if (m_gtol) feLog("    gap   : %15le %15le\n", gap, m_gtol);
    else        feLog("    gap   : %15le        ***\n", gap);
    double qctn = 0;
    if (fabs(normM1) > 1e-10) qctn = fabs((normM1 - normM0)/normM1);
    if (m_atol) feLog("    moment: %15le %15le\n", qctn, m_atol);
    else        feLog("    moment: %15le        ***\n", qctn);
    if (m_qtol) feLog("    angle : %15le %15le\n", qap, m_qtol);
    else        feLog("    angle : %15le        ***\n", qap);
    
    if (m_atol && ((pctn >= m_atol) || (qctn >= m_atol))) bconv = false;
    if (m_gtol && (gap >= m_gtol)) bconv = false;
    if (m_qtol && (qap >= m_qtol)) bconv = false;
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    if (!bconv)
    {
        // update multipliers
        m_L = Lm;
        m_U = Um;
    }
    
	if (m_bautopen)
	{
		// auto-penalty update (works only with gaptol and angtol)
		if (m_gtol && (gap > m_gtol)) {
			m_eps = fmax(gap / m_gtol, 100)*m_eps;
			feLog("    force_penalty :         %15le\n", m_eps);
		}
		if (m_qtol && (qap > m_qtol)) {
			m_ups = fmax(qap / m_qtol, 100)*m_ups;
			feLog("    moment_penalty :        %15le\n", m_ups);
		}
	}

    return bconv;
}

//-----------------------------------------------------------------------------
void FERigidLock::Update()
{
	if (m_laugon == FECore::LAGMULT_METHOD) return;

    vec3d ra, rb;
    vec3d za, zb;
    vec3d eat[3], eap[3], ea[3];
    vec3d ebt[3], ebp[3], eb[3];
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	const FETimeInfo& tp = GetTimeInfo();
	double alpha = tp.alphaf;
    
    ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*alpha + zap*(1-alpha);
	eat[0] = m_ea0[0]; RBa.GetRotation().RotateVector(eat[0]);
	eat[1] = m_ea0[1]; RBa.GetRotation().RotateVector(eat[1]);
	eat[2] = m_ea0[2]; RBa.GetRotation().RotateVector(eat[2]);
    eap[0] = m_ea0[0]; RBa.m_qp.RotateVector(eap[0]);
    eap[1] = m_ea0[1]; RBa.m_qp.RotateVector(eap[1]);
    eap[2] = m_ea0[2]; RBa.m_qp.RotateVector(eap[2]);
    ea[0] = eat[0]*alpha + eap[0]*(1-alpha);
    ea[1] = eat[1]*alpha + eap[1]*(1-alpha);
    ea[2] = eat[2]*alpha + eap[2]*(1-alpha);
    
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*alpha + zbp*(1-alpha);
	ebt[0] = m_eb0[0]; RBb.GetRotation().RotateVector(ebt[0]);
	ebt[1] = m_eb0[1]; RBb.GetRotation().RotateVector(ebt[1]);
	ebt[2] = m_eb0[2]; RBb.GetRotation().RotateVector(ebt[2]);
    ebp[0] = m_eb0[0]; RBb.m_qp.RotateVector(ebp[0]);
    ebp[1] = m_eb0[1]; RBb.m_qp.RotateVector(ebp[1]);
    ebp[2] = m_eb0[2]; RBb.m_qp.RotateVector(ebp[2]);
    eb[0] = ebt[0]*alpha + ebp[0]*(1-alpha);
    eb[1] = ebt[1]*alpha + ebp[1]*(1-alpha);
    eb[2] = ebt[2]*alpha + ebp[2]*(1-alpha);
    
    // incremental compound rotation of B w.r.t. A
    vec3d vth = ((ea[0] ^ eb[0]) + (ea[1] ^ eb[1]) + (ea[2] ^ eb[2]))/2;
    
    vec3d c = rb + zb - ra - za;
    m_F = m_L + c*m_eps;
    
    vec3d ksi = vth;
    m_M = m_U + ksi*m_ups;
    
}

//-----------------------------------------------------------------------------
void FERigidLock::Reset()
{
    m_F = vec3d(0,0,0);
    m_L = vec3d(0,0,0);
    m_M = vec3d(0,0,0);
    m_U = vec3d(0,0,0);
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    m_qa0 = m_q0 - RBa.m_r0;
    m_qb0 = m_q0 - RBb.m_r0;
    
    m_ea0[0] = m_e0[0]; m_ea0[1] = m_e0[1]; m_ea0[2] = m_e0[2];
    m_eb0[0] = m_e0[0]; m_eb0[1] = m_e0[1]; m_eb0[2] = m_e0[2];
}

//-----------------------------------------------------------------------------
vec3d FERigidLock::RelativeTranslation(const bool global)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;
    
    // body A
    vec3d ra = RBa.m_rt;
    vec3d za = m_qa0; RBa.GetRotation().RotateVector(za);
    
    // body B
    vec3d rb = RBb.m_rt;
    vec3d zb = m_qb0; RBb.GetRotation().RotateVector(zb);

    // relative translation in global coordinate system
    vec3d x = rb + zb - ra - za;
    
    if (global) return x;

    // evaluate local basis for body A
    vec3d ea[3];
    ea[0] = m_ea0[0]; RBa.GetRotation().RotateVector(ea[0]);
    ea[1] = m_ea0[1]; RBa.GetRotation().RotateVector(ea[1]);
    ea[2] = m_ea0[2]; RBa.GetRotation().RotateVector(ea[2]);

    // project relative translation onto local basis
    vec3d y(x*ea[0], x*ea[1], x*ea[2]);
    
    return y;
}

//-----------------------------------------------------------------------------
vec3d FERigidLock::RelativeRotation(const bool global)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;
    
    // get relative rotation
    quatd Q = RBb.GetRotation()*RBa.GetRotation().Inverse(); Q.MakeUnit();
    
    // relative rotation vector
    vec3d q = Q.GetRotationVector();
    
    if (global) return q;
    
    // evaluate local basis for body A
    vec3d ea[3];
    ea[0] = m_ea0[0]; RBa.GetRotation().RotateVector(ea[0]);
    ea[1] = m_ea0[1]; RBa.GetRotation().RotateVector(ea[1]);
    ea[2] = m_ea0[2]; RBa.GetRotation().RotateVector(ea[2]);

    // project relative rotation onto local basis
    vec3d y(q*ea[0], q*ea[1], q*ea[2]);
    
    return y;
}

void FERigidLock::Update(const std::vector<double>& ui)
{
	if (m_laugon == FECore::LAGMULT_METHOD)
	{
		m_F.x += ui[m_LM[0]];
		m_F.y += ui[m_LM[1]];
		m_F.z += ui[m_LM[2]];

		m_U1.x += ui[m_LM[3]];
		m_U1.y += ui[m_LM[4]];
		m_U1.z += ui[m_LM[5]];

		m_U2.x += ui[m_LM[6]];
		m_U2.y += ui[m_LM[7]];
		m_U2.z += ui[m_LM[8]];

		m_M = m_U1 + m_U2;
	}
}
