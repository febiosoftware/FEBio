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
#include "FERigidDamper.h"
#include "FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEMaterial.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidDamper, FERigidConnector);
	ADD_PARAMETER(m_c   , "c"          );
	ADD_PARAMETER(m_a0  , "insertion_a");
	ADD_PARAMETER(m_b0  , "insertion_b");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidDamper::FERigidDamper(FEModel* pfem) : FERigidConnector(pfem)
{
    m_nID = m_ncount++;
    m_c = 1.0;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidDamper::Init()
{
    // reset force
    m_F = vec3d(0,0,0);
    
	// base class first
	if (FERigidConnector::Init() == false) return false;

    // set spring insertions relative to rigid body center of mass
    m_qa0 = m_a0 - m_rbA->m_r0;
    m_qb0 = m_b0 - m_rbB->m_r0;
    
    m_at = m_a0;
    m_bt = m_b0;
    
    return true;
}

//-----------------------------------------------------------------------------
void FERigidDamper::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
	ar & m_qa0 & m_qb0;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidDamper::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<double> fa(6);
    vector<double> fb(6);
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	double alpha = tp.alphaf;
    
    // body A
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
    vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
    vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
    vec3d va = vat*alpha + vap*(1-alpha);
    
    // body b
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
    vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
    vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
    vec3d vb = vbt*alpha + vbp*(1-alpha);
    
    m_F = (vb - va)*m_c;
    
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
void FERigidDamper::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	double alpha = tp.alphaf;
	double beta  = tp.beta;
	double gamma = tp.gamma;
    
    // get time increment
    double dt = tp.timeIncrement;
    
    vector<int> LM(12);
	FEElementMatrix ke; ke.resize(12, 12);
    ke.zero();
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    mat3dd I(1);
    
    // body A
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d za = zat*alpha + zap*(1-alpha);
    mat3d zahat; zahat.skew(za);
    mat3d zathat; zathat.skew(zat);
    vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
    vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
    vec3d va = vat*alpha + vap*(1-alpha);
	quatd qai = RBa.GetRotation()*RBa.m_qp.Inverse(); qai.MakeUnit();
    vec3d cai = qai.GetVector()*(2*tan(qai.GetAngle()/2));
    mat3d Ta = I + skew(cai)/2 + dyad(cai)/4;
    
    // body b
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d zb = zbt*alpha + zbp*(1-alpha);
    mat3d zbhat; zbhat.skew(zb);
    mat3d zbthat; zbthat.skew(zbt);
    vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
    vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
    vec3d vb = vbt*alpha + vbp*(1-alpha);
	quatd qbi = RBb.GetRotation()*RBb.m_qp.Inverse(); qbi.MakeUnit();
    vec3d cbi = qbi.GetVector()*(2*tan(qbi.GetAngle()/2));
    mat3d Tb = I + skew(cbi)/2 + dyad(cbi)/4;
    
    m_F = (vb - va)*m_c;
    
    mat3ds A = I*(gamma/beta/dt);
    mat3d Ba = zathat*Ta.transpose()*(gamma/beta/dt) + skew(RBa.m_wt)*zathat;
    mat3d Bb = zbthat*Tb.transpose()*(gamma/beta/dt) + skew(RBb.m_wt)*zbthat;
    
    mat3d K;
    
    // (1,1)
    K = A*(alpha*m_c);
    ke[0][0] = K[0][0]; ke[0][1] = K[0][1]; ke[0][2] = K[0][2];
    ke[1][0] = K[1][0]; ke[1][1] = K[1][1]; ke[1][2] = K[1][2];
    ke[2][0] = K[2][0]; ke[2][1] = K[2][1]; ke[2][2] = K[2][2];
    
    // (1,2)
    K = Ba*(-alpha*m_c);
    ke[0][3] = K[0][0]; ke[0][4] = K[0][1]; ke[0][5] = K[0][2];
    ke[1][3] = K[1][0]; ke[1][4] = K[1][1]; ke[1][5] = K[1][2];
    ke[2][3] = K[2][0]; ke[2][4] = K[2][1]; ke[2][5] = K[2][2];
    
    // (1,3)
    K = A*(-alpha*m_c);
    ke[0][6] = K[0][0]; ke[0][7] = K[0][1]; ke[0][8] = K[0][2];
    ke[1][6] = K[1][0]; ke[1][7] = K[1][1]; ke[1][8] = K[1][2];
    ke[2][6] = K[2][0]; ke[2][7] = K[2][1]; ke[2][8] = K[2][2];
    
    // (1,4)
    K = Bb*(alpha*m_c);
    ke[0][9] = K[0][0]; ke[0][10] = K[0][1]; ke[0][11] = K[0][2];
    ke[1][9] = K[1][0]; ke[1][10] = K[1][1]; ke[1][11] = K[1][2];
    ke[2][9] = K[2][0]; ke[2][10] = K[2][1]; ke[2][11] = K[2][2];
    
    // (2,1)
    K = zahat*A*(alpha*m_c);
    ke[3][0] = K[0][0]; ke[3][1] = K[0][1]; ke[3][2] = K[0][2];
    ke[4][0] = K[1][0]; ke[4][1] = K[1][1]; ke[4][2] = K[1][2];
    ke[5][0] = K[2][0]; ke[5][1] = K[2][1]; ke[5][2] = K[2][2];
    
    // (2,2)
    K = (zahat*Ba)*(-alpha*m_c);
    ke[3][3] = K[0][0]; ke[3][4] = K[0][1]; ke[3][5] = K[0][2];
    ke[4][3] = K[1][0]; ke[4][4] = K[1][1]; ke[4][5] = K[1][2];
    ke[5][3] = K[2][0]; ke[5][4] = K[2][1]; ke[5][5] = K[2][2];
    
    // (2,3)
    K = zahat*A*(-alpha*m_c);
    ke[3][6] = K[0][0]; ke[3][7] = K[0][1]; ke[3][8] = K[0][2];
    ke[4][6] = K[1][0]; ke[4][7] = K[1][1]; ke[4][8] = K[1][2];
    ke[5][6] = K[2][0]; ke[5][7] = K[2][1]; ke[5][8] = K[2][2];
    
    // (2,4)
    K = (zahat*Bb)*(alpha*m_c);
    ke[3][9] = K[0][0]; ke[3][10] = K[0][1]; ke[3][11] = K[0][2];
    ke[4][9] = K[1][0]; ke[4][10] = K[1][1]; ke[4][11] = K[1][2];
    ke[5][9] = K[2][0]; ke[5][10] = K[2][1]; ke[5][11] = K[2][2];
    
    
    // (3,1)
    K = A*(-alpha*m_c);
    ke[6][0] = K[0][0]; ke[6][1] = K[0][1]; ke[6][2] = K[0][2];
    ke[7][0] = K[1][0]; ke[7][1] = K[1][1]; ke[7][2] = K[1][2];
    ke[8][0] = K[2][0]; ke[8][1] = K[2][1]; ke[8][2] = K[2][2];
    
    // (3,2)
    K = Ba*(alpha*m_c);
    ke[6][3] = K[0][0]; ke[6][4] = K[0][1]; ke[6][5] = K[0][2];
    ke[7][3] = K[1][0]; ke[7][4] = K[1][1]; ke[7][5] = K[1][2];
    ke[8][3] = K[2][0]; ke[8][4] = K[2][1]; ke[8][5] = K[2][2];
    
    // (3,3)
    K = A*(alpha*m_c);
    ke[6][6] = K[0][0]; ke[6][7] = K[0][1]; ke[6][8] = K[0][2];
    ke[7][6] = K[1][0]; ke[7][7] = K[1][1]; ke[7][8] = K[1][2];
    ke[8][6] = K[2][0]; ke[8][7] = K[2][1]; ke[8][8] = K[2][2];
    
    // (3,4)
    K = Bb*(-alpha*m_c);
    ke[6][9] = K[0][0]; ke[6][10] = K[0][1]; ke[6][11] = K[0][2];
    ke[7][9] = K[1][0]; ke[7][10] = K[1][1]; ke[7][11] = K[1][2];
    ke[8][9] = K[2][0]; ke[8][10] = K[2][1]; ke[8][11] = K[2][2];
    
    
    // (4,1)
    K = zbhat*A*(-alpha*m_c);
    ke[9 ][0] = K[0][0]; ke[ 9][1] = K[0][1]; ke[ 9][2] = K[0][2];
    ke[10][0] = K[1][0]; ke[10][1] = K[1][1]; ke[10][2] = K[1][2];
    ke[11][0] = K[2][0]; ke[11][1] = K[2][1]; ke[11][2] = K[2][2];
    
    // (4,2)
    K = (zbhat*Ba)*(alpha*m_c);
    ke[9 ][3] = K[0][0]; ke[ 9][4] = K[0][1]; ke[ 9][5] = K[0][2];
    ke[10][3] = K[1][0]; ke[10][4] = K[1][1]; ke[10][5] = K[1][2];
    ke[11][3] = K[2][0]; ke[11][4] = K[2][1]; ke[11][5] = K[2][2];
    
    // (4,3)
    K = zbhat*A*(alpha*m_c);
    ke[9 ][6] = K[0][0]; ke[ 9][7] = K[0][1]; ke[ 9][8] = K[0][2];
    ke[10][6] = K[1][0]; ke[10][7] = K[1][1]; ke[10][8] = K[1][2];
    ke[11][6] = K[2][0]; ke[11][7] = K[2][1]; ke[11][8] = K[2][2];
    
    // (4,4)
    K = (zbhat*Bb)*(-alpha*m_c);
    ke[9 ][9] = K[0][0]; ke[ 9][10] = K[0][1]; ke[ 9][11] = K[0][2];
    ke[10][9] = K[1][0]; ke[10][10] = K[1][1]; ke[10][11] = K[1][2];
    ke[11][9] = K[2][0]; ke[11][10] = K[2][1]; ke[11][11] = K[2][2];
    
    for (int j=0; j<6; ++j)
    {
        LM[j  ] = RBa.m_LM[j];
        LM[j+6] = RBb.m_LM[j];
    }
    
	ke.SetIndices(LM);
	LS.Assemble(ke);
}

//-----------------------------------------------------------------------------
bool FERigidDamper::Augment(int naug, const FETimeInfo& tp)
{
    return true;
}

//-----------------------------------------------------------------------------
void FERigidDamper::Update()
{
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    const FETimeInfo& tp = GetTimeInfo();
    double alpha = tp.alphaf;

    // body A
	vec3d zat = m_qa0; RBa.GetRotation().RotateVector(zat);
    vec3d zap = m_qa0; RBa.m_qp.RotateVector(zap);
    vec3d vat = RBa.m_vt + (RBa.m_wt ^ zat);
    vec3d vap = RBa.m_vp + (RBa.m_wp ^ zap);
    vec3d va = vat*alpha + vap*(1-alpha);
    
    // body b
	vec3d zbt = m_qb0; RBb.GetRotation().RotateVector(zbt);
    vec3d zbp = m_qb0; RBb.m_qp.RotateVector(zbp);
    vec3d vbt = RBb.m_vt + (RBb.m_wt ^ zbt);
    vec3d vbp = RBb.m_vp + (RBb.m_wp ^ zbp);
    vec3d vb = vbt*alpha + vbp*(1-alpha);
    
    m_at = RBa.m_rt + zat;
    m_bt = RBb.m_rt + zbt;
    
    m_F = (vb - va)*m_c;
}

//-----------------------------------------------------------------------------
void FERigidDamper::Reset()
{
    m_F = vec3d(0,0,0);
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    m_qa0 = m_a0 - RBa.m_r0;
    m_qb0 = m_b0 - RBb.m_r0;
}

//-----------------------------------------------------------------------------
vec3d FERigidDamper::RelativeTranslation(const bool global)
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

    return x;
}

//-----------------------------------------------------------------------------
vec3d FERigidDamper::RelativeRotation(const bool global)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;
    
    // get relative rotation
    quatd Q = RBb.GetRotation()*RBa.GetRotation().Inverse(); Q.MakeUnit();
    
    // relative rotation vector
    vec3d q = Q.GetRotationVector();
    
    return q;
}
