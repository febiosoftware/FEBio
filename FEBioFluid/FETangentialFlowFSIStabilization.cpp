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
#include "FETangentialFlowFSIStabilization.h"
#include "FEFluidMaterial.h"
#include "FEBioFluid.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FETangentialFlowFSIStabilization, FESurfaceLoad)
    ADD_PARAMETER(m_beta, "beta");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FETangentialFlowFSIStabilization::FETangentialFlowFSIStabilization(FEModel* pfem) : FESurfaceLoad(pfem), m_dofU(pfem), m_dofW(pfem)
{
    m_beta = 1.0;
    
    // get the degrees of freedom
    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofU.Clear();
        m_dofU.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::DISPLACEMENT));

        m_dofW.Clear();
        m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));

        m_dof.Clear();
        m_dof.AddDofs(m_dofU);
        m_dof.AddDofs(m_dofW);
    }
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETangentialFlowFSIStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FETangentialFlowFSIStabilization::Init()
{

    if (FESurfaceLoad::Init() == false) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
void FETangentialFlowFSIStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofU & m_dofW;
}

//-----------------------------------------------------------------------------
vec3d FETangentialFlowFSIStabilization::FluidVelocity(FESurfaceMaterialPoint& mp, double alpha)
{
    FESurfaceElement& el = *mp.SurfaceElement();
    vec3d vt[FEElement::MAX_NODES];
    int neln = el.Nodes();
    for (int j = 0; j<neln; ++j) {
        FENode& node = m_psurf->Node(el.m_lnode[j]);
        vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*alpha + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1 - alpha);
    }
    return el.eval(vt, mp.m_index);
}

//-----------------------------------------------------------------------------
void FETangentialFlowFSIStabilization::LoadVector(FEGlobalVector& R)
{
    const FETimeInfo& tp = GetTimeInfo();

    m_psurf->LoadVector(R, m_dofW, false, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {

        FESurfaceElement& el = *mp.SurfaceElement();
        
        // get the density
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
        double rho = fluid->ReferentialDensity();
        
        // tangent vectors
        vec3d rt[FEElement::MAX_NODES];
        m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);
        vec3d dxr = el.eval_deriv1(rt, mp.m_index);
        vec3d dxs = el.eval_deriv2(rt, mp.m_index);
        
        // normal and area element
        vec3d n = dxr ^ dxs;
        double da = n.unit();

        // fluid velocity
        vec3d v = FluidVelocity(mp, tp.alphaf);

        // tangential traction = -beta*density*|tangential velocity|*(tangential velocity)
        mat3dd I(1.0);
        vec3d vtau = (I - dyad(n))*v;
        double vmag = vtau.norm();

        // force vector (change sign for inflow vs outflow)
        vec3d f = vtau*(-m_beta*rho*vmag*da);

        double H = dof_a.shape;
        fa[0] = H * f.x;
        fa[1] = H * f.y;
        fa[2] = H * f.z;
    });
}

//-----------------------------------------------------------------------------
void FETangentialFlowFSIStabilization::StiffnessMatrix(FELinearSystem& LS)
{
    const FETimeInfo& tp = GetTimeInfo();

    m_psurf->LoadStiffness(LS, m_dof, m_dof, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
    
        FESurfaceElement& el = *mp.SurfaceElement();

        // fluid velocity
        vec3d v = FluidVelocity(mp, tp.alphaf);

        // get the density
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
        double rho = fluid->ReferentialDensity();
        
        // tangent vectors
        vec3d rt[FEElement::MAX_NODES];
        m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);
        vec3d dxr = el.eval_deriv1(rt, mp.m_index);
        vec3d dxs = el.eval_deriv2(rt, mp.m_index);
        
        vec3d n = dxr ^ dxs;
        double da = n.unit();
        
        mat3dd I(1.0);
        vec3d vtau = (I - dyad(n))*v;
        double vmag = vtau.unit();
        mat3d K = (I - dyad(n) + dyad(vtau))*(-m_beta*rho*vmag*da);
        // force vector (change sign for inflow vs outflow)
        vec3d ttt = vtau*(-m_beta*rho*vmag);

        // shape functions and derivatives
        double H_i  = dof_a.shape;

        double H_j  = dof_b.shape;
        double Gr_j = dof_b.shape_deriv_r;
        double Gs_j = dof_b.shape_deriv_s;

        // calculate stiffness component
        mat3d Kww = K*(H_i*H_j*tp.alphaf);
        vec3d g = (dxr*Gs_j - dxs*Gr_j)*(H_i*tp.alphaf);
        mat3d Kwu; Kwu.skew(g);
        Kwu = (ttt & n)*Kwu;
        Kab.zero();

        // dw/du
        Kab.sub(3, 0, Kwu);

        // dw/dw
        Kab.sub(3, 3, Kww);
    });
}
