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
#include "FERigidAngularDamper.h"
#include "FERigidBody.h"
#include "FECore/log.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEMaterial.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidAngularDamper, FERigidConnector);
	ADD_PARAMETER(m_c, "c");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidAngularDamper::FERigidAngularDamper(FEModel* pfem) : FERigidConnector(pfem)
{
    m_nID = m_ncount++;
}

//-----------------------------------------------------------------------------
//! TODO: This function is called twice: once in the Init and once in the Solve
//!       phase. Is that necessary?
bool FERigidAngularDamper::Init()
{
	// base class first
	if (FERigidConnector::Init() == false) return false;
    
    // reset force
    m_F = vec3d(0,0,0);
    
    return true;
}

//-----------------------------------------------------------------------------
void FERigidAngularDamper::Serialize(DumpStream& ar)
{
	FERigidConnector::Serialize(ar);
	ar & m_c;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidAngularDamper::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<double> fa(6);
    vector<double> fb(6);
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    double alpha = tp.alphaf;
    
    // body A
    vec3d wa = RBa.m_wt*alpha + RBa.m_wp*(1-alpha);
    
    // body b
    vec3d wb = RBb.m_wt*alpha + RBb.m_wp*(1-alpha);
    
    m_M = (wb - wa)*m_c;
    
    fa[0] = 0;
    fa[1] = 0;
    fa[2] = 0;
    
    fa[3] = m_M.x;
    fa[4] = m_M.y;
    fa[5] = m_M.z;
    
    fb[0] = 0;
    fb[1] = 0;
    fb[2] = 0;
    
    fb[3] = -m_M.x;
    fb[4] = -m_M.y;
    fb[5] = -m_M.z;
    
    for (int i=0; i<6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fa[i];
    for (int i=0; i<6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fb[i];
    
    RBa.m_Fr -= vec3d(fa[0],fa[1],fa[2]);
    RBa.m_Mr -= vec3d(fa[3],fa[4],fa[5]);
    RBb.m_Fr -= vec3d(fb[0],fb[1],fb[2]);
    RBb.m_Mr -= vec3d(fb[3],fb[4],fb[5]);
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FERigidAngularDamper::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    double alpha = tp.alphaf;
    double beta  = tp.beta;
    double gamma = tp.gamma;
    
    // get time increment
    double dt = tp.timeIncrement;
    
    int j;
    
    vector<int> LM(12);
	FEElementMatrix ke;
	ke.resize(12, 12);
    ke.zero();
    
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

    mat3dd I(1);
    
    // body A
    vec3d wa = RBa.m_wt*alpha + RBa.m_wp*(1-alpha);
	quatd qai = RBa.GetRotation()*RBa.m_qp.Inverse(); qai.MakeUnit();
    vec3d cai = qai.GetVector()*(2*tan(qai.GetAngle()/2));
    mat3d Ta = I + skew(cai)/2 + dyad(cai)/4;
    
    // body b
    vec3d wb = RBb.m_wt*alpha + RBb.m_wp*(1-alpha);
	quatd qbi = RBb.GetRotation()*RBb.m_qp.Inverse(); qbi.MakeUnit();
    vec3d cbi = qbi.GetVector()*(2*tan(qbi.GetAngle()/2));
    mat3d Tb = I + skew(cbi)/2 + dyad(cbi)/4;
    
    m_M = (wb - wa)*m_c;
    
    mat3d Ba = Ta.transpose()*(gamma/beta/dt);
    mat3d Bb = Tb.transpose()*(gamma/beta/dt);
    
    mat3d K;
    K.zero();
    
    // (2,2)
    K = Ba*(alpha*m_c);
    ke[3][3] = K[0][0]; ke[3][4] = K[0][1]; ke[3][5] = K[0][2];
    ke[4][3] = K[1][0]; ke[4][4] = K[1][1]; ke[4][5] = K[1][2];
    ke[5][3] = K[2][0]; ke[5][4] = K[2][1]; ke[5][5] = K[2][2];
    
    // (2,4)
    K = Bb*(-alpha*m_c);
    ke[3][9] = K[0][0]; ke[3][10] = K[0][1]; ke[3][11] = K[0][2];
    ke[4][9] = K[1][0]; ke[4][10] = K[1][1]; ke[4][11] = K[1][2];
    ke[5][9] = K[2][0]; ke[5][10] = K[2][1]; ke[5][11] = K[2][2];
    
    
    // (4,2)
    K = Ba*(-alpha*m_c);
    ke[9 ][3] = K[0][0]; ke[ 9][4] = K[0][1]; ke[ 9][5] = K[0][2];
    ke[10][3] = K[1][0]; ke[10][4] = K[1][1]; ke[10][5] = K[1][2];
    ke[11][3] = K[2][0]; ke[11][4] = K[2][1]; ke[11][5] = K[2][2];
    
    // (4,4)
    K = Bb*(alpha*m_c);
    ke[9 ][9] = K[0][0]; ke[ 9][10] = K[0][1]; ke[ 9][11] = K[0][2];
    ke[10][9] = K[1][0]; ke[10][10] = K[1][1]; ke[10][11] = K[1][2];
    ke[11][9] = K[2][0]; ke[11][10] = K[2][1]; ke[11][11] = K[2][2];
    
    for (j=0; j<6; ++j)
    {
        LM[j  ] = RBa.m_LM[j];
        LM[j+6] = RBb.m_LM[j];
    }
    
	ke.SetIndices(LM);
	LS.Assemble(ke);
}

//-----------------------------------------------------------------------------
bool FERigidAngularDamper::Augment(int naug, const FETimeInfo& tp)
{
    return true;
}

//-----------------------------------------------------------------------------
void FERigidAngularDamper::Update()
{
	FERigidBody& RBa = *m_rbA;
	FERigidBody& RBb = *m_rbB;

	const FETimeInfo& tp = GetTimeInfo();
	double alpha = tp.alphaf;
    
    // body A
    vec3d wa = RBa.m_wt*alpha + RBa.m_wp*(1-alpha);
    
    // body b
    vec3d wb = RBb.m_wt*alpha + RBb.m_wp*(1-alpha);
    
    m_M = (wb - wa)*m_c;
}

//-----------------------------------------------------------------------------
void FERigidAngularDamper::Reset()
{
    m_F = vec3d(0,0,0);
    m_M = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
vec3d FERigidAngularDamper::RelativeTranslation(const bool global)
{
    return vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
vec3d FERigidAngularDamper::RelativeRotation(const bool global)
{
    FERigidBody& RBa = *m_rbA;
    FERigidBody& RBb = *m_rbB;
    
    // get relative rotation
    quatd Q = RBb.GetRotation()*RBa.GetRotation().Inverse(); Q.MakeUnit();
    
    // relative rotation vector
    vec3d q = Q.GetRotationVector();
    
    return q;
}
