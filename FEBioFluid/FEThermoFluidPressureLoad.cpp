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
#include "FEThermoFluidPressureLoad.h"
#include "FEBioThermoFluid.h"
#include "FEThermoFluidSolver.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEThermoFluidPressureLoad, FESurfaceLoad)
    ADD_PARAMETER(m_p, "pressure")->setUnits("P")->setLongName("fluid pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEThermoFluidPressureLoad::FEThermoFluidPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_p = 0;
}

//-----------------------------------------------------------------------------
bool FEThermoFluidPressureLoad::Init()
{
    
    m_dofEF = m_dof.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION));
    m_dofT = m_dof.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE));

    m_Pr = GetGlobalConstant("P");
    
    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! serialization
void FEThermoFluidPressureLoad::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofT & m_dofEF;
}

//-----------------------------------------------------------------------------
bool FEThermoFluidPressureLoad::FluidPressure(FEElement& el, double& p, double& dpJ, double& dpT)
{
    FEMaterial* pm = GetFEModel()->GetMaterial(el.GetMatID());
    FEThermoFluid* fluid = pm->ExtractProperty<FEThermoFluid>();
    if (fluid == nullptr) return false;
    int nint = el.GaussPoints();
    p = dpJ = dpT = 0;
    for (int n=0; n<nint; ++n) {
        FEMaterialPoint* mp = el.GetMaterialPoint(n);
        FEFluidMaterialPoint* fp = mp->ExtractData<FEFluidMaterialPoint>();
        p += fluid->GetElastic()->Pressure(*mp);
        dpJ += fluid->GetElastic()->Tangent_Strain(*mp);
        dpT += fluid->GetElastic()->Tangent_Temperature(*mp);
    }
    p /= nint;
    dpJ /= nint;
    dpT /= nint;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEThermoFluidPressureLoad::FluidPressure2ndDerivs(FEElement& el, double& dpJJ, double& dpJT, double& dpTT)
{
    FEMaterial* pm = GetFEModel()->GetMaterial(el.GetMatID());
    FEThermoFluid* fluid = pm->ExtractProperty<FEThermoFluid>();
    if (fluid == nullptr) return false;
    int nint = el.GaussPoints();
    dpJJ = dpJT = dpTT = 0;
    for (int n=0; n<nint; ++n) {
        FEMaterialPoint* mp = el.GetMaterialPoint(n);
        FEFluidMaterialPoint* fp = mp->ExtractData<FEFluidMaterialPoint>();
        dpJJ += fluid->GetElastic()->Tangent_Strain_Strain(*mp);
        dpJT += fluid->GetElastic()->Tangent_Strain_Temperature(*mp);
        dpTT += fluid->GetElastic()->Tangent_Temperature_Temperature(*mp);
    }
    dpJJ /= nint;
    dpJT /= nint;
    dpTT /= nint;
    
    return true;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FEThermoFluidPressureLoad::LoadVector(FEGlobalVector& R)
{
    double alpha = GetFEModel()->GetTime().alphaf;
    double dt = GetFEModel()->GetTime().timeIncrement;
    FESurface& surf = GetSurface();

    surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {
        double p0 = m_p(pt);
        // get the underlying thermofluid element
        FESurfaceElement& el = *pt.SurfaceElement();
        FEElement* pe = el.m_elem[0].pe;
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEThermoFluid* fluid = pm->ExtractProperty<FEThermoFluid>();
        if (fluid == nullptr) return;
        double p, dpJ, dpT;
        if (!FluidPressure(*pe,p,dpJ,dpT)) return;
        double da = (pt.dxr ^ pt.dxs).norm();
        double Na  = dof_a.shape;
        
        val[0] = Na*dpJ*(p - p0)*da;
        val[1] = Na*dpT*(p - p0)*da;
    });

}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FEThermoFluidPressureLoad::StiffnessMatrix(FELinearSystem& LS)
{
    double alpha = GetFEModel()->GetTime().alphaf;
    double dt = GetFEModel()->GetTime().timeIncrement;
    FESurface& surf = GetSurface();
    surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& kab) {

        double p0 = m_p(mp);
        // get the underlying thermofluid element
        FESurfaceElement& el = *mp.SurfaceElement();
        FEElement* pe = el.m_elem[0].pe;
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEThermoFluid* fluid = pm->ExtractProperty<FEThermoFluid>();
        if (fluid == nullptr) return;
        double p, dpJ, dpT;
        if (!FluidPressure(*pe,p,dpJ,dpT)) return;
        double dpJJ, dpJT, dpTT;
        if (!FluidPressure2ndDerivs(*pe,dpJJ,dpJT,dpTT)) return;

        double da = (mp.dxr ^ mp.dxs).norm();
        double Na  = dof_a.shape;
        double Nb  = dof_b.shape;

        kab(0, 0) = -(Na*Nb*(dpJJ*(p-p0)+pow(dpJ,2))*da);
        kab(0, 1) = -(Na*Nb*(dpJT*(p-p0)+dpJ*dpT)*da);
        kab(1, 0) = kab(0, 1);
        kab(1, 1) = -(Na*Nb*(dpTT*(p-p0)+pow(dpT,2))*da);
    });
}
