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
#include "FESoluteNaturalFlux.h"
#include "FEMultiphasic.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESoluteNaturalFlux, FESurfaceLoad)
    ADD_PARAMETER(m_bshellb, "shell_bottom");
    ADD_PARAMETER(m_isol   , "solute_id");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESoluteNaturalFlux::FESoluteNaturalFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofC(pfem), m_dofU(pfem), m_dofP(pfem)
{
    m_bshellb = false;
    m_isol = 0;
}
    
//-----------------------------------------------------------------------------
//! allocate storage
void FESoluteNaturalFlux::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! serialization
void FESoluteNaturalFlux::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);

    if (ar.IsShallow() == false)
    {
        ar & m_dofC & m_dofU & m_dofP;
    }
}

//-----------------------------------------------------------------------------
bool FESoluteNaturalFlux::Init()
{
    if (m_isol == -1) return false;

    // set up the dof lists
    FEModel* fem = GetFEModel();
    m_dofC.Clear();
    m_dofU.Clear();
    m_dofP.Clear();
    if (m_bshellb == false)
    {
        m_dofC.AddDof(fem->GetDOFIndex("concentration", m_isol - 1));

        m_dofU.AddDof(fem->GetDOFIndex("x"));
        m_dofU.AddDof(fem->GetDOFIndex("y"));
        m_dofU.AddDof(fem->GetDOFIndex("z"));
        
        m_dofP.AddDof(fem->GetDOFIndex("p"));
    }
    else
    {
        m_dofC.AddDof(fem->GetDOFIndex("shell concentration", m_isol - 1));

        m_dofU.AddDof(fem->GetDOFIndex("sx"));
        m_dofU.AddDof(fem->GetDOFIndex("sy"));
        m_dofU.AddDof(fem->GetDOFIndex("sz"));

        m_dofP.AddDof(fem->GetDOFIndex("q"));
    }
    m_dof.AddDofs(m_dofU);
    m_dof.AddDofs(m_dofP);
    m_dof.AddDofs(m_dofC);
    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FESoluteNaturalFlux::LoadVector(FEGlobalVector& R)
{
    double dt = CurrentTimeIncrement();

    m_psurf->SetShellBottom(m_bshellb);

    FESoluteNaturalFlux* flux = this;
    m_psurf->LoadVector(R, m_dofC, false, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, std::vector<double>& fa) {

        // get surface element
        FESurfaceElement& el = *mp.SurfaceElement();
        // get underlying solid element
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        // get the local solute id
        FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
        if (psi == nullptr) {
            fa[0] = 0;
            return;
        }
        int sid = psi->FindLocalSoluteID(flux->m_isol);
        
        // get element-averaged fluid flux and actual solute concentration
        vec3d w(0,0,0);
        double c = 0;
        int nint = pe->GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(pt.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
            w += pb.m_w;
            c += ps.m_ca[sid];
        }
        w /= nint;
        c /= nint;

        // evaluate desired natural solute flux
        vec3d dxt = mp.dxr ^ mp.dxs;
        double jn = c*(w*dxt);
        if (flux->m_bshellb) jn = -jn;

        // molar flow rate
        double f = jn* dt;

        double H_i = dof_a.shape;
        fa[0] = H_i * f;
    });
}

//-----------------------------------------------------------------------------
void FESoluteNaturalFlux::StiffnessMatrix(FELinearSystem& LS)
{
    // time increment
    double dt = CurrentTimeIncrement();

    m_psurf->SetShellBottom(m_bshellb);
    
    // evaluate the stiffness contribution
    FESoluteNaturalFlux* flux = this;
    m_psurf->LoadStiffness(LS, m_dofC, m_dofU, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {

        // get surface element
        FESurfaceElement& el = *mp.SurfaceElement();
        // get underlying solid element
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        // get the local solute id
        FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
        if (psi == nullptr) return;
        int sid = psi->FindLocalSoluteID(flux->m_isol);
        
        // get element-averaged fluid flux and actual solute concentration
        vec3d w(0,0,0);
        double c = 0;
        int nint = pe->GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(pt.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
            w += pb.m_w;
            c += ps.m_ca[sid];
        }
        w /= nint;
        c /= nint;
        
        // shape functions and derivatives
        double H_i  = dof_a.shape;
        double Gr_j = dof_b.shape_deriv_r;
        double Gs_j = dof_b.shape_deriv_s;

        // calculate surface normal
        vec3d dxt = mp.dxr ^ mp.dxs;
        vec3d nu = dxt.normalized();
        double jn = c*(w*nu);
        if (flux->m_bshellb) jn = -jn;

        // calculate stiffness component
        vec3d t1 = nu*jn;
        vec3d t2 = mp.dxs*Gr_j - mp.dxr*Gs_j;
        vec3d kab = (t1 ^ t2)*(H_i)*dt;

        Kab[0][0] = kab.x;
        Kab[0][1] = kab.y;
        Kab[0][2] = kab.z;
    });
}
