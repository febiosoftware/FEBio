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
#include "FERigidFollowerMoment.h"
#include "FERigidBody.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEMaterial.h"
#include "FECore/FELoadCurve.h"
#include "FEMechModel.h"
#include "FERigidMaterial.h"
#include <FECore/FELinearSystem.h>

//=============================================================================
BEGIN_FECORE_CLASS(FERigidFollowerMoment, FERigidLoad);
    ADD_PARAMETER(m_rid      , "rb"       )->setEnums("$(rigid_materials)");
    ADD_PARAMETER(m_m        , "moment"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidFollowerMoment::FERigidFollowerMoment(FEModel* pfem) : FERigidLoad(pfem)
{
    m_rid = -1;
    m_m = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
//! do some sanity checks
bool FERigidFollowerMoment::Init()
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
void FERigidFollowerMoment::Serialize(DumpStream& ar)
{
    FEModelLoad::Serialize(ar);
    ar & m_rid;
    ar & m_m;
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidFollowerMoment::LoadVector(FEGlobalVector& R)
{
    FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidBody& body = *fem.GetRigidBody(m_rid);

    const FETimeInfo& tp = fem.GetTime();
    double alpha = tp.alphaf;
    
    // calculate the moment value
    vec3d mt = body.GetRotation()*m_m;
    vec3d mp = body.GetPreviousRotation()*m_m;
    vec3d m = mt*alpha + mp*(1-alpha);

    // apply force and moment to body
    int n;
    n = body.m_LM[3]; if (n >= 0) R[n] += m.x;
    n = body.m_LM[4]; if (n >= 0) R[n] += m.y;
    n = body.m_LM[5]; if (n >= 0) R[n] += m.z;
    
    body.m_Mr += m;
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
//! TODO: Only the stiffness contribution in the were the axial forces are applied
//!       to the center of mass has been implemented.
void FERigidFollowerMoment::StiffnessMatrix(FELinearSystem& LS)
{
    FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
    FERigidBody& body = *fem.GetRigidBody(m_rid);

    const FETimeInfo& tp = fem.GetTime();
    double alpha = tp.alphaf;
    
    // calculate the moment value
    vec3d mt = body.GetRotation()*m_m;

    // build the stiffness matrix components
    mat3da mthat(mt);
    mat3d Kqq = mthat*(-alpha);

    // put it all together
    FEElementMatrix K; K.resize(3, 3); K.zero();
    K.sub(0,0, Kqq);

    // get the equation numbers
    vector<int> lm(3);
    lm[ 0] = body.m_LM[3];
    lm[ 1] = body.m_LM[4];
    lm[ 2] = body.m_LM[5];

    // assemble into global matrix
    K.SetIndices(lm);
    LS.Assemble(K);
}
