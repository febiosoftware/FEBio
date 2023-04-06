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
#include "FEFluidNaturalHeatFlux.h"
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FEGlobalVector.h>
#include <FECore/log.h>
#include <FECore/LinearSolver.h>
#include <FECore/FEModel.h>
#include "FEBioThermoFluid.h"
#include "FEFluidMaterialPoint.h"
#include "FEThermoFluidMaterialPoint.h"

//-----------------------------------------------------------------------------
//! constructor
FEFluidNaturalHeatFlux::FEFluidNaturalHeatFlux(FEModel* pfem) : FESurfaceLoad(pfem)
{
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidNaturalHeatFlux::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidNaturalHeatFlux::Init()
{
    // get the degrees of freedom
    m_dof.Clear();
    m_dof.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE));
    if (m_dof.IsEmpty()) return false;

    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the prescribed normal velocity
void FEFluidNaturalHeatFlux::LoadVector(FEGlobalVector& R)
{
    m_psurf->LoadVector(R, m_dof, true, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, std::vector<double>& fa) {
        
        const FETimeInfo& tp = GetTimeInfo();
        
        // get surface element
        FESurfaceElement& el = *mp.SurfaceElement();
        // get underlying solid element
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        vec3d dxt = mp.dxr ^ mp.dxs;

        // get element-averaged heat flux
        vec3d q(0,0,0);
        int nint = pe->GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEFluidMaterialPoint* pf = pt.ExtractData<FEFluidMaterialPoint>();
            if (pf == nullptr) break;
            FEThermoFluidMaterialPoint* ps = pt.ExtractData<FEThermoFluidMaterialPoint>();
            if (ps == nullptr) break;
//            if (ps->m_q*pf->m_vft > 0) q += ps->m_q;
            q += ps->m_q;
        }
        q /= nint;
        
        double qn = dxt*q;
        double H = dof_a.shape;
        
        fa[0] = qn*H;
    });
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidNaturalHeatFlux::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
}
