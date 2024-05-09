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
#include "FEFluidMovingFrameLoad.h"
#include "FEFluidMaterialPoint.h"
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEFluidMovingFrameLoad, FEBodyForce)
    ADD_PROPERTY(m_ksi[0], "qx");
    ADD_PROPERTY(m_ksi[1], "qy");
    ADD_PROPERTY(m_ksi[2], "qz");
    ADD_PROPERTY(m_ct[0], "cx");
    ADD_PROPERTY(m_ct[1], "cy");
    ADD_PROPERTY(m_ct[2], "cz");
END_FECORE_CLASS();

FEFluidMovingFrameLoad::FEFluidMovingFrameLoad(FEModel* fem) : FEBodyForce(fem)
{
    m_omega = m_alpha = m_c = m_cdot = m_cddot = vec3d(0,0,0);
    m_q = quatd(vec3d(0,0,0));
}

bool FEFluidMovingFrameLoad::Init()
{
    for (int k=0; k<3; ++k) {
        if (!m_ksi[k]->Init()) return false;
        if (!m_ct[k]->Init()) return false;
    }
    
    return FEBodyForce::Init();
}

void FEFluidMovingFrameLoad::PrepStep()
{
    const FETimeInfo& tp = GetTimeInfo();
    double dt = tp.timeIncrement;
    double t = tp.currentTime - (1-tp.alphaf)*dt;

    vec3d ksi(m_ksi[0]->value(t), m_ksi[1]->value(t), m_ksi[2]->value(t));
    m_q = quatd(ksi);
    m_q.MakeUnit();
    m_omega = vec3d(m_ksi[0]->derive(t), m_ksi[1]->derive(t), m_ksi[2]->derive(t));
    m_alpha = vec3d(m_ksi[0]->deriv2(t), m_ksi[1]->deriv2(t), m_ksi[2]->deriv2(t));
    m_c = vec3d(m_ct[0]->value(t), m_ct[1]->value(t), m_ct[2]->value(t));
    m_cdot = vec3d(m_ct[0]->derive(t), m_ct[1]->derive(t), m_ct[2]->derive(t));
    m_cddot = vec3d(m_ct[0]->deriv2(t), m_ct[1]->deriv2(t), m_ct[2]->deriv2(t));
}

vec3d FEFluidMovingFrameLoad::force(FEMaterialPoint& pt)
{
    FEFluidMaterialPoint& fp = *pt.ExtractData<FEFluidMaterialPoint>();
    
    vec3d x = pt.m_r0;
    m_q.RotateVector(x);
    vec3d v = fp.m_vft;
    m_q.RotateVector(v);

    // FEBio negates the body force, so we use the negative of the value in the notes
    vec3d f = (m_cddot + (m_alpha ^ x) + (m_omega ^ (m_omega ^ x)) + (m_omega ^ v)*2);
    m_q.Inverse().RotateVector(f);

    return f;
}

mat3d FEFluidMovingFrameLoad::stiffness(FEMaterialPoint& pt)
{
    const FETimeInfo& tp = GetTimeInfo();

    mat3d What; What.skew(m_omega);

    // for fluid analyses the position is not variable
    mat3d K = m_q.Inverse().RotationMatrix()*What*m_q.RotationMatrix()*(2*tp.alphaf);

    return K;
}
