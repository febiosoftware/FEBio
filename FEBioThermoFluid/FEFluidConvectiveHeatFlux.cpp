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
#include "FEFluidConvectiveHeatFlux.h"
#include "FEBioThermoFluid.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FEModel.h>
#include <FECore/matrix.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidConvectiveHeatFlux, FESurfaceLoad)
    ADD_PARAMETER(m_h    , "h")->setLongName("Constant heat transfer coefficient");
    ADD_PARAMETER(m_Tinf , "Tinf")->setLongName("Sink temperature");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidConvectiveHeatFlux::FEFluidConvectiveHeatFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dof(pfem)
{
//    m_h = 0;
//    m_Tinf = 0;
    m_dofT = pfem->GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE), 0);
    m_dof.AddDof(m_dofT);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidConvectiveHeatFlux::Init()
{
    // get the degrees of freedom
    m_dof.Clear();
    m_dof.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE));
    if (m_dof.IsEmpty()) return false;

    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the prescribed convective heat flux
void FEFluidConvectiveHeatFlux::LoadVector(FEGlobalVector& R)
{
    FESurface& surf = GetSurface();
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();

    // evaluate the integral
    surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {
        
        double T[FEElement::MAX_NODES];
        double Tavg = 0;

        // get surface element
        FESurfaceElement& el = *pt.SurfaceElement();
        for (int i=0; i<el.Nodes(); ++i) {
            FENode& node = mesh.Node(el.m_node[i]);
            Tavg += node.get(m_dofT);
        }
        Tavg /= el.Nodes();
        
        // evaluate heat flux at this material point
        double h = m_h(pt);
        double Tinf = m_Tinf(pt);

        double q = h*(Tavg - Tinf);

        double J = (pt.dxr ^ pt.dxs).norm();

        double N = dof_a.shape;

        val[0] = -N*q*J;
    });
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness matrix for the prescribed convective heat flux
void FEFluidConvectiveHeatFlux::StiffnessMatrix(FELinearSystem& LS)
{
    FESurface& surf = GetSurface();
    const FETimeInfo& tp = GetTimeInfo();

    surf.LoadStiffness(LS, m_dof, m_dof, [=](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
        
        // shape functions and derivatives
        double H_i  = dof_a.shape;
        
        double H_j  = dof_b.shape;
        double Gr_j = dof_b.shape_deriv_r;
        double Gs_j = dof_b.shape_deriv_s;
        
        double h = m_h(pt);
        double J = (pt.dxr ^ pt.dxs).norm();
        matrix kabhT(1,1);
        kabhT(0,0) = H_i*H_j*h*J*tp.alphaf;
        
        Kab = kabhT;
    });
}
