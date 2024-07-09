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
#include "FETemperatureStabilization.h"
#include "FEFluid.h"
#include "FEBioThermoFluid.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FETemperatureStabilization, FESurfaceLoad)
    ADD_PARAMETER(m_beta, "beta");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FETemperatureStabilization::FETemperatureStabilization(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_beta = 1.0;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETemperatureStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FETemperatureStabilization::Init()
{
    FEModel& fem = *GetFEModel();
    
    // determine the nr of concentration equations
    m_dofW.Clear();
    m_dofW.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::RELATIVE_FLUID_VELOCITY));
    m_dofT = fem.GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE), 0);
    m_dof.Clear();
    m_dof.AddDofs(m_dofW);
    m_dof.AddDof(m_dofT);

    FESurface* ps = &GetSurface();

    m_nnlist.Create(fem.GetMesh());

    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! calculate residual
void FETemperatureStabilization::LoadVector(FEGlobalVector& R)
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
        vec3d r0[FEElement::MAX_NODES];
        m_psurf->GetReferenceNodalCoordinates(el, r0);
        vec3d dxr = el.eval_deriv1(r0, mp.m_index);
        vec3d dxs = el.eval_deriv2(r0, mp.m_index);
        
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
void FETemperatureStabilization::StiffnessMatrix(FELinearSystem& LS)
{
    const FETimeInfo& tp = GetTimeInfo();

    m_psurf->LoadStiffness(LS, m_dofW, m_dofW, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
    
        FESurfaceElement& el = *mp.SurfaceElement();

        // get the density
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
        double rho = fluid->ReferentialDensity();

        // fluid velocity
        vec3d v = FluidVelocity(mp, tp.alphaf);

        // tangent vectors
        vec3d r0[FEElement::MAX_NODES];
        m_psurf->GetReferenceNodalCoordinates(el, r0);
        vec3d dxr = el.eval_deriv1(r0, mp.m_index);
        vec3d dxs = el.eval_deriv2(r0, mp.m_index);
        
        vec3d n = dxr ^ dxs;
        double da = n.unit();
        
        mat3dd I(1.0);
        vec3d vtau = (I - dyad(n))*v;
        double vmag = vtau.unit();
        mat3d K = (I - dyad(n) + dyad(vtau))*(-m_beta*rho*vmag*da);

        // shape functions and derivatives
        double H_i  = dof_a.shape;

        double H_j  = dof_b.shape;

        // calculate stiffness component
        mat3d Kww = K*(H_i*H_j*tp.alphaf);
        Kab.zero();

        // dw/dw
        Kab.sub(0, 0, Kww);
    });
}
//-----------------------------------------------------------------------------
//! serialization
void FETemperatureStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
	ar & m_dofW;
	ar & m_dofT;
    m_nnlist.Create(GetFEModel()->GetMesh());
}
