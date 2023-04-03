/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include "FETangentialFlowPFStabilization.h"
#include "FEBioPolarFluid.h"
#include "FEPolarFluid.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FETangentialFlowPFStabilization, FESurfaceLoad)
    ADD_PARAMETER(m_beta, "beta");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FETangentialFlowPFStabilization::FETangentialFlowPFStabilization(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem), m_dofG(pfem)
{
    m_beta = 1.0;
    
    // get the degrees of freedom
    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofW.Clear();
        m_dofW.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::RELATIVE_FLUID_VELOCITY));

        m_dofG.Clear();
        m_dofG.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_ANGULAR_VELOCITY));
        
        m_dof.Clear();
        m_dof.AddDofs(m_dofW);
        m_dof.Clear();
        m_dof.AddDofs(m_dofW);
        m_dof.AddDofs(m_dofG);
    }
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETangentialFlowPFStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FETangentialFlowPFStabilization::Init()
{
    
    if (FESurfaceLoad::Init() == false) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
void FETangentialFlowPFStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofW & m_dofG;
}

//-----------------------------------------------------------------------------
vec3d FETangentialFlowPFStabilization::FluidVelocity(FESurfaceMaterialPoint& mp, double alpha)
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
vec3d FETangentialFlowPFStabilization::FluidAngularVelocity(FESurfaceMaterialPoint& mp, double alpha)
{
    FESurfaceElement& el = *mp.SurfaceElement();
    vec3d gt[FEElement::MAX_NODES];
    int neln = el.Nodes();
    for (int j = 0; j<neln; ++j) {
        FENode& node = m_psurf->Node(el.m_lnode[j]);
        gt[j] = node.get_vec3d(m_dofG[0], m_dofG[1], m_dofG[2])*alpha + node.get_vec3d_prev(m_dofG[0], m_dofG[1], m_dofG[2])*(1 - alpha);
    }
    return el.eval(gt, mp.m_index);
}

//-----------------------------------------------------------------------------
void FETangentialFlowPFStabilization::LoadVector(FEGlobalVector& R)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    
    m_psurf->LoadVector(R, m_dof, false, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {
        
        FESurfaceElement& el = *mp.SurfaceElement();
        
        // get the density
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEPolarFluid* pfluid = pm->ExtractProperty<FEPolarFluid>();
        FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
        double rho = fluid->ReferentialDensity();
        double kg = (pfluid) ? pfluid->m_kg : 0;
        
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
        
        // fluid angular velocity
        vec3d g = FluidAngularVelocity(mp, tp.alphaf);
        
        // tangential traction = -beta*density*|tangential velocity|*(tangential velocity)
        mat3dd I(1.0);
        vec3d vtau = (I - dyad(n))*v;
        double vmag = vtau.norm();
        vec3d gtau = (I - dyad(n))*g;
        double gmag = gtau.norm();

        // force vector (change sign for inflow vs outflow)
        vec3d fv = vtau*(-m_beta*rho*vmag*da);
        vec3d fg = gtau*(-m_beta*rho*pow(kg,3)*gmag*da);

        double H = dof_a.shape;
        fa[0] = H * fv.x;
        fa[1] = H * fv.y;
        fa[2] = H * fv.z;
        fa[3] = H * fg.x;
        fa[4] = H * fg.y;
        fa[5] = H * fg.z;
    });
}

//-----------------------------------------------------------------------------
void FETangentialFlowPFStabilization::StiffnessMatrix(FELinearSystem& LS)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    
    m_psurf->LoadStiffness(LS, m_dof, m_dof, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
        
        FESurfaceElement& el = *mp.SurfaceElement();
        
        // fluid velocity
        vec3d v = FluidVelocity(mp, tp.alphaf);
        
        // fluid angular velocity
        vec3d g = FluidAngularVelocity(mp, tp.alphaf);
        
        // get the density
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEPolarFluid* pfluid = pm->ExtractProperty<FEPolarFluid>();
        FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
        double rho = fluid->ReferentialDensity();
        double kg = (pfluid) ? pfluid->m_kg : 0;

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
        mat3d Kv = (I - dyad(n) + dyad(vtau))*(-m_beta*rho*vmag*da);
        
        vec3d gtau = (I - dyad(n))*g;
        double gmag = gtau.unit();
        mat3d Kg = (I - dyad(n) + dyad(gtau))*(-m_beta*rho*pow(kg,3)*gmag*da);
        
        // shape functions and derivatives
        double H_i  = dof_a.shape;
        
        double H_j  = dof_b.shape;

        // calculate stiffness component
        mat3d Kww = Kv*(H_i*H_j*tp.alphaf);
        mat3d Kgg = Kg*(H_i*H_j*tp.alphaf);
        Kab.zero();
        
        // dw/dw
        Kab.sub(0, 0, Kww);
        // dg/dg
        Kab.sub(3, 3, Kgg);
    });
}
