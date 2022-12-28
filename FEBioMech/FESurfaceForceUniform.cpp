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
#include "FESurfaceForceUniform.h"
#include "FEBioMech.h"
#include <FECore/FEFacetSet.h>

//=============================================================================
BEGIN_FECORE_CLASS(FESurfaceForceUniform, FESurfaceLoad)
    ADD_PARAMETER(m_scale   , "scale");
    ADD_PARAMETER(m_force   , "force");
    ADD_PARAMETER(m_bshellb , "shell_bottom");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESurfaceForceUniform::FESurfaceForceUniform(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_scale = 1.0;
    m_force = m_traction = vec3d(0, 0, 0);
    m_bshellb = false;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FESurfaceForceUniform::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
// initialization
bool FESurfaceForceUniform::Init()
{
    // get the degrees of freedom
    m_dof.Clear();
    if (m_bshellb == false)
    {
        m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
    }
    else
    {
        m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
    }
    if (m_dof.IsEmpty()) return false;
    
    if (FESurfaceLoad::Init() == false) return false;
    
    // evaluate uniform traction based on applied force
    FESurface& surf = GetSurface();
    surf.SetShellBottom(m_bshellb);
    double area = 0;
    for (int i=0; i<surf.Elements(); ++i) {
        FESurfaceElement& el = surf.Element(i);
        area += surf.FaceArea(el);
    }
    m_traction = m_force/area;
    
    return true;
}

//-----------------------------------------------------------------------------
void FESurfaceForceUniform::LoadVector(FEGlobalVector& R)
{
    FESurface& surf = GetSurface();
    surf.SetShellBottom(m_bshellb);
    
    // evaluate the integral
    FESurfaceForceUniform* load = this;
    surf.LoadVector(R, m_dof, true, [=](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {
        
        // evaluate traction at this material point
        vec3d t = m_traction*m_scale;
        if (load->m_bshellb) t = -t;
        
        double J = (pt.dxr ^ pt.dxs).norm();
        
        double H_u = dof_a.shape;
        
        val[0] = H_u * t.x*J;
        val[1] = H_u * t.y*J;
        val[2] = H_u * t.z*J;
    });
}

//-----------------------------------------------------------------------------
void FESurfaceForceUniform::StiffnessMatrix(FELinearSystem& LS)
{
    // Nothing to do here.
    // TODO: I think if the linear flag is false, I do need to evaluate a stiffness.
}
