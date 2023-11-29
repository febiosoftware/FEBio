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

//For now only allows for one solute at a time, jn=-d0*gradc+c*vf
//Stiffness does not include effect from different solutes



#include "stdafx.h"
#include "FEFluidSolutesNaturalFlux.h"
#include "FEFluidMaterial.h"
#include "FEFluidSolutes.h"
#include "FEBioFluidSolutes.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include "FEFluidSolutesDomain3D.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEFluidSolutesNaturalFlux, FESurfaceLoad)
    ADD_PARAMETER(m_isol   , "solute_id")->setEnums("$(solutes)");
    ADD_PARAMETER(m_bup    , "update");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEFluidSolutesNaturalFlux::FEFluidSolutesNaturalFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem), m_dofEF(pfem), m_dofC(pfem)
{
    m_isol = -1;
    m_bup = false;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidSolutesNaturalFlux::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidSolutesNaturalFlux::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    
    if (ar.IsShallow() == false)
    {
        ar & m_dofC & m_dofW & m_dofEF;
    }
}

//-----------------------------------------------------------------------------
bool FEFluidSolutesNaturalFlux::Init()
{
    if (m_isol <= 0) return false;
    
    // set up the dof lists
    FEModel* fem = GetFEModel();
    m_dofC.Clear();
    m_dofW.Clear();
    m_dofEF.Clear();
    m_dofC.AddDof(fem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), m_isol - 1));
    m_dofW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY));
    m_dofEF.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION));

    m_dof.AddDofs(m_dofW);
    m_dof.AddDofs(m_dofC);

    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FEFluidSolutesNaturalFlux::Update()
{
    if (m_bup) {
        for (int is=0; is<m_psurf->Elements(); ++is)
        {
            // get surface element
            FESurfaceElement& el = m_psurf->Element(is);
            // get underlying solid element
            FESolidElement* pe = dynamic_cast<FESolidElement*>(el.m_elem[0]);
            if (pe == nullptr) break;
            // get element data
            int neln = pe->Nodes();
            int nint = pe->GaussPoints();
            
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            // get the local solute id
            FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
            if (psi == nullptr) break;
            int sid = psi->FindLocalSoluteID(m_isol);
            if (sid == -1) break;
            
            // identify nodes on the surface
            vector<bool> nsrf(neln,false);
            for (int j=0; j<neln; ++j) {
                for (int k=0; k < el.Nodes(); ++k) {
                    if (el.m_node[k] == pe->m_node[j]) nsrf[j] = true;
                }
            }
            
            // get average effective concentration of nodes not on surface
            double cavg = 0;
            int m = 0;
            for (int i=0; i<neln; ++i) {
                if (!nsrf[i]) {
                    int n = pe->m_node[i];
                    FENode& node = GetMesh().Node(n);
                    int dof = m_dofC[m_isol-1];
                    if (dof != -1) {
                        cavg += node.get(dof);
                        ++m;
                    }
                }
            }
            // assign this average value to surface nodes as initial guess
            if (m) {
                cavg /= m;
                for (int i=0; i<neln; ++i) {
                    if (nsrf[i]) {
                        int n = pe->m_node[i];
                        FENode& node = GetMesh().Node(n);
                        int dof = m_dofC[m_isol-1];
                        if (dof != -1) node.set(dof, cavg);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesNaturalFlux::LoadVector(FEGlobalVector& R)
{
    FEFluidSolutesNaturalFlux* flux = this;
    m_psurf->LoadVector(R, m_dofC, true, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, std::vector<double>& fa) {
        
        const FETimeInfo& tp = GetTimeInfo();
        
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
        FEFluidSolutes* pmat = dynamic_cast<FEFluidSolutes*>(pm);
        int sid = psi->FindLocalSoluteID(flux->m_isol);
        vec3d dxt = mp.dxr ^ mp.dxs;
        
        // get element-averaged diffusive flux
        vec3d jd(0,0,0);
        int nint = pe->GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            jd += pmat->SoluteDiffusiveFlux(pt, m_isol-1);
        }
        jd /= nint;
        
        // evaluate desired natural solute flux = normal convective flux * area
        double jn = jd*dxt*tp.alphaf;
        
        
        double H_i = dof_a.shape;
        fa[0] = H_i * jn;
    });
}
/*{
    FEFluidSolutesNaturalFlux* flux = this;
    m_psurf->LoadVector(R, m_dofC, true, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, std::vector<double>& fa) {
        
        const FETimeInfo& tp = GetTimeInfo();
        
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
        vec3d dxt = mp.dxr ^ mp.dxs;

        // get element-averaged partition coefficient
        double kappa = 0;
        int nint = pe->GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEFluidSolutesMaterialPoint& ps = *(pt.ExtractData<FEFluidSolutesMaterialPoint>());
            kappa += ps.m_k[sid];
        }
        kappa /= nint;

        // evaluate average effective solute concentration and fluid velocity at nodes
        vec3d w(0,0,0);
        double ce = 0, cp = 0;
        for (int i=0; i<el.Nodes(); ++i) {
            FENode& node = m_psurf->Node(el.m_lnode[i]);
            w += node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*tp.alphaf + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1-tp.alphaf);
            ce += node.get(m_dofC[m_isol-1])*tp.alphaf + node.get_prev(m_dofC[m_isol-1])*(1-tp.alphaf);
            cp += node.get_prev(m_dofC[m_isol-1]);
        }
        w /= el.Nodes();
        ce /= el.Nodes();
        cp /= el.Nodes();

        // evaluate desired natural solute flux = normal convective flux * area
        double wn = w*dxt;
        double jn = kappa*ce*wn;

        // molar flow rate (if negative, use previous value of concentration)
        double f = (wn >= 0) ? jn : kappa*cp*wn;
//        double f = (jn > 0) ? jn : 0;

        double H_i = dof_a.shape;
        fa[0] = H_i * f;
    });
}*/

//-----------------------------------------------------------------------------
void FEFluidSolutesNaturalFlux::StiffnessMatrix(FELinearSystem& LS)
{}
/*{
    FEFluidSolutesNaturalFlux* flux = this;
    m_psurf->LoadStiffness(LS, m_dofC, m_dof, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
        
        const FETimeInfo& tp = GetTimeInfo();
        
        // get surface element
        FESurfaceElement& el = *mp.SurfaceElement();
        // get underlying solid element
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        // get the local solute id
        FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
        if (psi == nullptr) return;
        int sid = psi->FindLocalSoluteID(flux->m_isol);
        
        // get element-averaged solute partition coefficient
        double kappa = 0, d0 = 0;
        int nint = pe->GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEFluidSolutesMaterialPoint& ps = *(pt.ExtractData<FEFluidSolutesMaterialPoint>());
            kappa += ps.m_k[sid];
            d0 += psi->GetSolute(sid)->m_pDiff->Free_Diffusivity(pt);
        }
        kappa /= nint;
        d0 /= nint;

        // evaluate average effective solute concentration and fluid velocity
        vec3d w(0,0,0);
        double ce = 0;
        for (int i=0; i<el.Nodes(); ++i) {
            FENode& node = m_psurf->Node(el.m_lnode[i]);
            w += node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*tp.alphaf + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1-tp.alphaf);
            ce += node.get(m_dofC[m_isol-1])*tp.alphaf + node.get_prev(m_dofC[m_isol-1])*(1-tp.alphaf);
        }
        w /= el.Nodes();
        ce /= el.Nodes();

        // shape functions and derivatives
        double H_i  = dof_a.shape;
        double H_j  = dof_b.shape;
        
        // calculate surface normal
        vec3d dxt = mp.dxr ^ mp.dxs;
        
        // calculate stiffness component
        vec3d kcv = dxt*(H_i*H_j*kappa*ce*tp.alphaf);
        double kcc = (dxt*w)*(H_i*H_j*kappa)*tp.alphaf;
        
        Kab[0][0] = -kcv.x;
        Kab[0][1] = -kcv.y;
        Kab[0][2] = -kcv.z;
        Kab[0][3] = -kcc;
    });
}*/
