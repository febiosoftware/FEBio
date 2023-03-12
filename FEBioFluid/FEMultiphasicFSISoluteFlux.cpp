/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEMultiphasicFSISoluteFlux.h"
#include "FEBioMultiphasicFSI.h"
#include <FECore/log.h>
#include <FECore/FEFacetSet.h>
#include "FECore/FEAnalysis.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEMultiphasicFSISoluteFlux, FESurfaceLoad)
ADD_PARAMETER(m_flux   , "flux");
ADD_PARAMETER(m_isol   , "sol");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEMultiphasicFSISoluteFlux::FEMultiphasicFSISoluteFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofC(pfem), m_dofU(pfem)
{
    m_flux = 1.0;
    m_isol = -1;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEMultiphasicFSISoluteFlux::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
    m_flux.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the prescribed normal component of velocity
void FEMultiphasicFSISoluteFlux::LoadVector(FEGlobalVector& R)
{
    FEMultiphasicFSISoluteFlux* flux = this;
    m_psurf->LoadVector(R, m_dofC, true, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {
        
        double wr = flux->m_flux(mp);
        
        vec3d dxt = mp.dxr ^ mp.dxs;
        
        // volumetric flow rate
        double f = dxt.norm()*wr;
        
        double H_i = dof_a.shape;
        fa[0] = H_i * f;
    });
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSISoluteFlux::StiffnessMatrix(FELinearSystem& LS)
{
    const FETimeInfo& tp = GetTimeInfo();

    // evaluate the stiffness contribution
    FEMultiphasicFSISoluteFlux* flux = this;
    m_psurf->LoadStiffness(LS, m_dofC, m_dofU, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
        
        // shape functions and derivatives
        double H_i  = dof_a.shape;
        double Gr_j = dof_b.shape_deriv_r;
        double Gs_j = dof_b.shape_deriv_s;
        
        double wr = flux->m_flux(mp);
        
        // calculate surface normal
        vec3d dxt = mp.dxr ^ mp.dxs;
        
        double alpha = tp.alphaf;
        
        // calculate stiffness component
        vec3d t1 = dxt / dxt.norm()*wr;
        vec3d t2 = - mp.dxs*Gr_j + mp.dxr*Gs_j;
        mat3d A; A.skew(t2);
        vec3d kab = -A*t1*(H_i)*alpha;
        Kab.zero();
        Kab[0][0] -= kab.x;
        Kab[0][1] -= kab.y;
        Kab[0][2] -= kab.z;
    });
}

//-----------------------------------------------------------------------------
//! initialize
bool FEMultiphasicFSISoluteFlux::Init()
{
    if (m_isol == -1) return false;
    
    // set up the dof lists
    FEModel* fem = GetFEModel();
    m_dofC.Clear();
    m_dofU.Clear();
    m_dofC.AddDof(GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_CONCENTRATION), m_isol-1));
    m_dofU.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::DISPLACEMENT));
    m_dof.Clear();
    m_dof.AddDofs(m_dofU);
    m_dof.AddDofs(m_dofC);
    
    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! serialization
void FEMultiphasicFSISoluteFlux::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    
    if (ar.IsShallow() == false)
    {
        ar & m_dofC & m_dofU;
    }
}
