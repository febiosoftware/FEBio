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
#include "FERigidFollowerForce.h"
#include "FERigidBody.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEMaterial.h"
#include "FECore/FELoadCurve.h"
#include "FEMechModel.h"
#include "FERigidMaterial.h"
#include <FECore/FELinearSystem.h>

//=============================================================================
BEGIN_FECORE_CLASS(FERigidFollowerForce, FERigidLoad);
    ADD_PARAMETER(m_rid      , "rb"       )->setEnums("$(rigid_materials)");
    ADD_PARAMETER(m_X        , "insertion");
    ADD_PARAMETER(m_f        , "force"    );
    ADD_PARAMETER(m_brelative, "relative" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidFollowerForce::FERigidFollowerForce(FEModel* pfem) : FERigidLoad(pfem)
{
    m_rid = -1;
    m_X = m_f = vec3d(0,0,0);
    m_brelative = false;
}

//-----------------------------------------------------------------------------
//! do some sanity checks
bool FERigidFollowerForce::Init()
{
    // At this point the rigid ID's are still associated with the materials.
    // We want to associate them with the rigid objects instead.
    FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_rid-1));
    if (pm == 0) return false;
    m_rid = pm->GetRigidBodyID(); if (m_rid < 0) return false;

    // all is well in the world
    return true;
}

//-----------------------------------------------------------------------------
void FERigidFollowerForce::Serialize(DumpStream& ar)
{
    FEModelLoad::Serialize(ar);
    ar & m_rid;
    ar & m_X & m_f;
    ar & m_brelative;
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidFollowerForce::LoadVector(FEGlobalVector& R)
{
    FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidBody& body = *fem.GetRigidBody(m_rid);
    const FETimeInfo& tp = fem.GetTime();

    double alpha = tp.alphaf;
    
    // get the attachment position in global coordinates for body A
    vec3d Z = (m_brelative ? m_X : m_X - body.m_r0);
    vec3d zt = body.GetRotation()*Z;
    vec3d zp = body.GetPreviousRotation()*Z;
    vec3d z = zt*alpha + zp*(1-alpha);

    // calculate the force value
    vec3d ft = body.GetRotation()*m_f;
    vec3d fp = body.GetPreviousRotation()*m_f;
    vec3d f = ft*alpha + fp*(1-alpha);

    // get the moment about the center of mass
    vec3d m = z ^ f;

    // apply force and moment to body
    int n;
    n = body.m_LM[0]; if (n >= 0) R[n] += f.x;
    n = body.m_LM[1]; if (n >= 0) R[n] += f.y;
    n = body.m_LM[2]; if (n >= 0) R[n] += f.z;
    n = body.m_LM[3]; if (n >= 0) R[n] += m.x;
    n = body.m_LM[4]; if (n >= 0) R[n] += m.y;
    n = body.m_LM[5]; if (n >= 0) R[n] += m.z;
    
    body.m_Fr += f;
    body.m_Mr += m;
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
//! TODO: Only the stiffness contribution in the were the axial forces are applied
//!       to the center of mass has been implemented.
void FERigidFollowerForce::StiffnessMatrix(FELinearSystem& LS)
{
    FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidBody& body = *fem.GetRigidBody(m_rid);

    const FETimeInfo& tp = fem.GetTime();
    double alpha = tp.alphaf;
    
    // get the attachment position in global coordinates for body A
    vec3d Z = (m_brelative ? m_X : m_X - body.m_r0);
    vec3d zt = body.GetRotation()*Z;
    vec3d zp = body.GetPreviousRotation()*Z;
    vec3d z = zt*alpha + zp*(1-alpha);

    // calculate the force value
    vec3d ft = body.GetRotation()*m_f;
    vec3d fp = body.GetPreviousRotation()*m_f;
    vec3d f = ft*alpha + fp*(1-alpha);

    // build the stiffness matrix components
    mat3da fthat(ft);
    mat3da zthat(zt);
    mat3da fhat(f);
    mat3da zhat(z);
    mat3d Krq = fthat*(-alpha);
    mat3d Kqq = (fhat*zthat - zhat*fthat)*alpha;

    // put it all together
    FEElementMatrix K; K.resize(6, 6); K.zero();
    K.sub(0,3, Krq);
    K.sub(3,3, Kqq);

    // get the equation numbers
    vector<int> lm(6);
    lm[ 0] = body.m_LM[0];
    lm[ 1] = body.m_LM[1];
    lm[ 2] = body.m_LM[2];
    lm[ 3] = body.m_LM[3];
    lm[ 4] = body.m_LM[4];
    lm[ 5] = body.m_LM[5];

    // assemble into global matrix
    K.SetIndices(lm);
    LS.Assemble(K);
}
