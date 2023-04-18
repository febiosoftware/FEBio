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
#include "FEFluidSolutesFlux.h"
#include "FEBioFluidSolutes.h"
#include <FEBioMix/FESoluteInterface.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FEAnalysis.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidSolutesFlux, FESurfaceLoad)
ADD_PARAMETER(m_flux   , "flux")->setUnits("n/L^2.t")->setLongName("effective solute molar flux");
ADD_PARAMETER(m_isol   , "solute_id")->setEnums("$(solutes)");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidSolutesFlux::FEFluidSolutesFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofC(pfem)
{
    m_flux = 1.0;
    m_isol = -1;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidSolutesFlux::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
    m_flux.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the prescribed normal component of velocity
void FEFluidSolutesFlux::LoadVector(FEGlobalVector& R)
{
    FEFluidSolutesFlux* flux = this;
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
//! initialize
bool FEFluidSolutesFlux::Init()
{
    if (m_isol == -1) return false;
    
    // set up the dof lists
    FEModel* fem = GetFEModel();
    m_dofC.Clear();
    
    m_dofC.AddDof(fem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), m_isol-1));

    m_dof.Clear();
    m_dof.AddDofs(m_dofC);
    
    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidSolutesFlux::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    
    if (ar.IsShallow() == false)
    {
        ar & m_dofC;
    }
}
