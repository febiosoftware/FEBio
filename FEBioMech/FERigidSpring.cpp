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
#include "FERigidSpring.h"
#include "FERigidBody.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidSpring, FERigidConnector);
    ADD_PARAMETER(m_k   , "k"          );
    ADD_PARAMETER(m_a0  , "insertion_a");
    ADD_PARAMETER(m_b0  , "insertion_b");
    ADD_PARAMETER(m_L0  , "free_length");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidSpring::FERigidSpring(FEModel* pfem) : FERigidConnector(pfem)
{
    m_nID = m_ncount++;
    m_L0 = 0;
    m_k = 1.0;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidSpring::Init()
{
	// base class first
	if (FERigidConnector::Init() == false) return false;
    
    // reset force
    m_F = vec3d(0,0,0);
    
    // if not specified by use, get free length of spring
    if (m_L0 == 0) m_L0 = (m_b0 - m_a0).norm();

    // set spring insertions relative to rigid body center of mass
    m_qa0 = m_a0 - m_rbA->m_r0;
    m_qb0 = m_b0 - m_rbB->m_r0;

	m_at = m_a0;
	m_bt = m_b0;
    
    return true;
}

//-----------------------------------------------------------------------------
void FERigidSpring::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
    ar & m_qa0 & m_qb0;
    ar & m_L0 & m_k;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidSpring::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<double> fa(6);
    vector<double> fb(6);
    
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;

	double alpha = tp.alphaf;

    // body A
	quatd Qat = RBa.GetRotation();
	quatd Qap = RBa.m_qp;
    vec3d ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
	vec3d zat = Qat*m_qa0;
    vec3d zap = Qap*m_qa0;
    vec3d za = zat*alpha + zap*(1-alpha);
    
    // body b
	quatd Qbt = RBb.GetRotation();
	quatd Qbp = RBb.m_qp;
    vec3d rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
	vec3d zbt = Qbt*m_qb0;
    vec3d zbp = Qbp*m_qb0;
    vec3d zb = zbt*alpha + zbp*(1-alpha);
    
    vec3d c = rb + zb - ra - za;
    double L = c.norm();
    m_F = (L > 0) ? c*((1 - m_L0/L)*m_k) : c*m_k;
    
    fa[0] = m_F.x;
    fa[1] = m_F.y;
    fa[2] = m_F.z;
    
    fa[3] = za.y*m_F.z-za.z*m_F.y;
    fa[4] = za.z*m_F.x-za.x*m_F.z;
    fa[5] = za.x*m_F.y-za.y*m_F.x;
    
    fb[0] = -m_F.x;
    fb[1] = -m_F.y;
    fb[2] = -m_F.z;
    
    fb[3] = -zb.y*m_F.z+zb.z*m_F.y;
    fb[4] = -zb.z*m_F.x+zb.x*m_F.z;
    fb[5] = -zb.x*m_F.y+zb.y*m_F.x;
    
    for (int i=0; i<6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
    for (int i=0; i<6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];
    
    RBa.m_Fr -= vec3d(fa[0],fa[1],fa[2]);
    RBa.m_Mr -= vec3d(fa[3],fa[4],fa[5]);
    RBb.m_Fr -= vec3d(fb[0],fb[1],fb[2]);
    RBb.m_Mr -= vec3d(fb[3],fb[4],fb[5]);
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidSpring::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	double alpha = tp.alphaf;
    
    vector<int> LM(12);
    FEElementMatrix ke(12,12);
    ke.zero();
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    // body A
	quatd Qat = RBa.GetRotation();
	quatd Qap = RBa.m_qp;
	vec3d rat = RBa.m_rt;
	vec3d rap = RBa.m_rp;
    vec3d ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
    mat3d zahat; zahat.skew(za);
    mat3d zathat; zathat.skew(zat);
    
    // body b
	quatd Qbt = RBb.GetRotation();
	quatd Qbp = RBb.m_qp;
	vec3d rbt = RBb.m_rt;
	vec3d rbp = RBb.m_rp;
	vec3d rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
    mat3d zbhat; zbhat.skew(zb);
    mat3d zbthat; zbthat.skew(zbt);
    
    vec3d c = rb + zb - ra - za;
    double L = c.norm();
    vec3d n = c/L;
    m_F = (L > 0) ? c*((1 - m_L0/L)*m_k) : c*m_k;
    mat3ds P;
    mat3dd I(1);
    P = (L > 0) ? I*(1 - m_L0/L) + dyad(n)*(m_L0/L) : I;

	mat3da Fhat(m_F);
    
    mat3d K;
    
    // row 0
	ke.set(0, 0, P * (m_k));
	ke.set(0, 3, P * zathat * (-m_k));
	ke.set(0, 6, P * (-m_k));
	ke.set(0, 9, P * zbthat * (m_k));
    
    // row 1
	ke.set(3, 0, (zahat * P) * (m_k));
	ke.set(3, 3, (zahat * P * zathat) * (-m_k) + (Fhat * zathat) * (-1.0));
	ke.set(3, 6, (zahat * P) * (- m_k));
	ke.set(3, 9, (zahat * P * zbthat) * (m_k));
    
    // row 2
	ke.set(6, 0, P * (-m_k));
	ke.set(6, 3, P * zathat * (m_k));
	ke.set(6, 6, P * (m_k));
	ke.set(6, 9, P * zbthat * (-m_k));
    
    // row 3
	ke.set(9, 0, (zbhat * P)* (-m_k));
	ke.set(9, 3, (zbhat * P * zathat) * (m_k));
	ke.set(9, 6, (zbhat * P)* (m_k));
	ke.set(9, 9, (zbhat * P * zbthat) * (-m_k) + (Fhat * zbthat));

	if (alpha != 1) ke *= alpha;
    
    for (int j=0; j<6; ++j)
    {
        LM[j  ] = RBa.m_LM[j];
        LM[j+6] = RBb.m_LM[j];
    }
    
	ke.SetIndices(LM);
	LS.Assemble(ke);
}

//-----------------------------------------------------------------------------
bool FERigidSpring::Augment(int naug, const FETimeInfo& tp)
{
    return true;
}

//-----------------------------------------------------------------------------
void FERigidSpring::Update()
{
    vec3d ra, rb, c;
    vec3d za, zb;
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	const FETimeInfo& tp = GetTimeInfo();
	double alpha = tp.alphaf;

    ra = RBa.m_rt*alpha + RBa.m_rp*(1-alpha);
    rb = RBb.m_rt*alpha + RBb.m_rp*(1-alpha);
    
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    za = zat*alpha + zap*(1-alpha);
    
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    zb = zbt*alpha + zbp*(1-alpha);

	m_at = RBa.m_rt + zat;
	m_bt = RBb.m_rt + zbt;

    c = rb + zb - ra - za;
    double L = c.norm();
    m_F = (L > 0) ? c*((1 - m_L0/L)*m_k) : c*m_k;
}

//-----------------------------------------------------------------------------
void FERigidSpring::Reset()
{
    m_F = vec3d(0,0,0);
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    m_qa0 = m_a0 - RBa.m_r0;
    m_qb0 = m_b0 - RBb.m_r0;
}

//-----------------------------------------------------------------------------
vec3d FERigidSpring::RelativeTranslation(const bool global)
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
    
    vec3d y = vec3d(1,0,0)*x.norm();

    return y;
}

//-----------------------------------------------------------------------------
vec3d FERigidSpring::RelativeRotation(const bool global)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;
    
    // get relative rotation
    quatd Q = RBb.GetRotation()*RBa.GetRotation().Inverse(); Q.MakeUnit();
    
    // relative rotation vector
    vec3d q = Q.GetRotationVector();
    
    return q;
}
