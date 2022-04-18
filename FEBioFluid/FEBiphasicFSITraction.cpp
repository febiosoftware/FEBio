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
#include "FEBiphasicFSITraction.h"
#include "FECore/FEModel.h"
#include "FEFluid.h"
#include "FEBiphasicFSI.h"
#include "FEBioFSI.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEBiphasicFSITraction, FESurfaceLoad)
ADD_PARAMETER(m_bshellb , "shell_bottom");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEBiphasicFSITraction::FEBiphasicFSITraction(FEModel* pfem) : FESurfaceLoad(pfem), m_dofU(pfem), m_dofSU(pfem), m_dofW(pfem)
{
    m_bshellb = false;

    // get the degrees of freedom
    // TODO: Can this be done in Init, since  there is no error checking
    if (pfem)
    {
        m_dofU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT));
        m_dofSU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_DISPLACEMENT));
        m_dofW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY));
        m_dofEF = pfem->GetDOFIndex(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION), 0);

        m_dof.Clear();
        m_dof.AddDofs(m_dofU);
        m_dof.AddDofs(m_dofSU);
        m_dof.AddDofs(m_dofW);
        m_dof.AddDof(m_dofEF);
    }
}

//-----------------------------------------------------------------------------
//! initialize
bool FEBiphasicFSITraction::Init()
{
    FESurface& surf = GetSurface();
    surf.SetShellBottom(m_bshellb);
    surf.SetInterfaceStatus(true);
    if (FESurfaceLoad::Init() == false) return false;
    
    // get the list of fluid-FSI elements connected to this interface
    FEModel* fem = GetFEModel();
    int NF = surf.Elements();
    m_elem.resize(NF);
    m_s.resize(NF, 1);
    for (int j = 0; j<NF; ++j)
    {
        bool bself = false;
        FESurfaceElement& el = surf.Element(j);
        // extract the first of two elements on this interface
        m_elem[j] = el.m_elem[0];
        if (el.m_elem[1] == nullptr) bself = true;
        // get its material and check if FEBiphasicFSI
        FEMaterial* pm = fem->GetMaterial(m_elem[j]->GetMatID());
        FEBiphasicFSI* pfsi = dynamic_cast<FEBiphasicFSI*>(pm);
        if (pfsi) {
            double s = m_psurf->FacePointing(el, *m_elem[j]);
            m_s[j] = bself ? -s : s;
            if (m_s[j] == 0) return false;
        }
        else if (!bself) {
            // extract the second of two elements on this interface
            m_elem[j] = el.m_elem[1];
            pm = fem->GetMaterial(m_elem[j]->GetMatID());
            pfsi = dynamic_cast<FEBiphasicFSI*>(pm);
            if (pfsi == nullptr) return false;
            m_s[j] = m_psurf->FacePointing(el, *m_elem[j]);
            if (m_s[j] == 0) return false;
        }
        else
            return false;
    }
    
    // TODO: Deal with the case when the surface is a shell domain separating two FSI domains
    // that use different fluid bulk moduli
    
    return true;
}

//-----------------------------------------------------------------------------
double FEBiphasicFSITraction::GetFluidDilatation(FESurfaceMaterialPoint& mp, double alpha)
{
    double ef = 0;
    FESurfaceElement& el = *mp.SurfaceElement();
    double* H = el.H(mp.m_index);
    int neln = el.Nodes();
    for (int j = 0; j < neln; ++j) {
        FENode& node = m_psurf->Node(el.m_lnode[j]);
        double ej = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1.0 - alpha);
        ef += ej*H[j];
    }
    return ef;
}

//-----------------------------------------------------------------------------
mat3ds FEBiphasicFSITraction::GetFluidStress(FESurfaceMaterialPoint& pt)
{
    FEModel* fem = GetFEModel();
    FESurfaceElement& face = *pt.SurfaceElement();
    int iel = face.m_lid;
    
    // Get the fluid stress from the fluid-FSI element
    mat3ds sv(mat3dd(0));
    FEElement* pe = m_elem[iel];
    int nint = pe->GaussPoints();
    FEBiphasicFSI* pfsi = dynamic_cast<FEBiphasicFSI*>(fem->GetMaterial(pe->GetMatID()));
    for (int n = 0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
        sv += pfsi->Fluid()->GetViscous()->Stress(mp);
    }
    sv /= nint;
    return sv;
}

//-----------------------------------------------------------------------------
void FEBiphasicFSITraction::LoadVector(FEGlobalVector& R)
{
    const FETimeInfo& tp = GetTimeInfo();

    // If surface is bottom of shell, we should take shell displacement dofs (i.e. m_dofSU).
    FEDofList dof = m_bshellb ? m_dofSU : m_dofU;
    m_psurf->LoadVector(R, dof, false, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {
        
        // get the surface element
        FESurfaceElement& el = *mp.SurfaceElement();
        int iel = el.m_lid;
        FEBiphasicFSI* pfsi = dynamic_cast<FEBiphasicFSI*>(GetFEModel()->GetMaterial(m_elem[iel]->GetMatID()));

        // nodal coordinates
        vec3d rt[FEElement::MAX_NODES];
        m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);
        
        // evaluate covariant basis vectors at integration point
        vec3d gr = el.eval_deriv1(rt, mp.m_index)*m_s[iel];
        vec3d gs = el.eval_deriv2(rt, mp.m_index);
        vec3d gt = gr ^ gs;
        
        // Get the fluid stress at integration point
        // necessarily using the attached solid element
        mat3ds sv = GetFluidStress(mp);
        
        // fluid dilatation at integration point
        // only from surface element
        double ef = GetFluidDilatation(mp, tp.alphaf);
        double p = pfsi->Fluid()->Pressure(ef);

        // evaluate traction
        vec3d f = gt*p - sv*gt;
        
        double H = dof_a.shape;
        fa[0] = H * f.x;
        fa[1] = H * f.y;
        fa[2] = H * f.z;
    });
}

//-----------------------------------------------------------------------------
void FEBiphasicFSITraction::StiffnessMatrix(FELinearSystem& LS)
{
    const FETimeInfo& tp = GetTimeInfo();

    FEModel* fem = GetFEModel();
    FESurface* ps = &GetSurface();
    
    // build dof list
    // TODO: If surface is bottom of shell, we should take shell displacement dofs (i.e. m_dofSU).
    FEDofList dofs(fem);
    if (!m_bshellb) dofs.AddDofs(m_dofU); else dofs.AddDofs(m_dofSU);
    dofs.AddDofs(m_dofW);
    dofs.AddDof(m_dofEF);
    
    // evaluate stiffness
    m_psurf->LoadStiffness(LS, dofs, dofs, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
        
        FESurfaceElement& el = *mp.SurfaceElement();
        int iel = el.m_lid;
        int neln = el.Nodes();
        
        double dt = tp.timeIncrement;
        double alpha = tp.alphaf;
        double a = tp.gamma / (tp.beta*dt);
        
        vector<vec3d> gradN(neln);
        
        // nodal coordinates
        vec3d rt[FEElement::MAX_NODES];
        ps->GetNodalCoordinates(el, tp.alphaf, rt);
        
        // Get the fluid stress and its tangents from the fluid-FSI element
        mat3ds sv(mat3dd(0)), svJ(mat3dd(0));
        tens4ds cv; cv.zero();
        mat3d Ls; Ls.zero();
        mat3d Dw; Dw.zero();
        double phif=0.0;
        double phis=0.0;
        double J=0.0;
        vec3d gradphif = vec3d(0.0);
        vec3d gradJ = vec3d(0.0);
        FEElement* pe = m_elem[iel];
        int pint = pe->GaussPoints();
        FEBiphasicFSI* pfsi = dynamic_cast<FEBiphasicFSI*>(fem->GetMaterial(pe->GetMatID()));
        for (int n = 0; n<pint; ++n)
        {
            FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
            FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
            FEBiphasicFSIMaterialPoint& ft = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
            sv += pfsi->Fluid()->GetViscous()->Stress(mp);
            svJ += pfsi->Fluid()->GetViscous()->Tangent_Strain(mp);
            cv += pfsi->Fluid()->Tangent_RateOfDeformation(mp);
            Ls += ep.m_L;
            phif += pfsi->Porosity(mp);
            phis += pfsi->SolidVolumeFrac(mp);
            gradphif += pfsi->gradPorosity(mp);
            gradJ += ft.m_gradJ;
            Dw += ft.m_Lw.sym();
            J += ep.m_J;
        }
        sv /= pint;
        svJ /= pint;
        cv /= pint;
        Ls /= pint;
        phif /= pint;
        phis /=pint;
        gradphif /=pint;
        gradJ /=pint;
        Dw /=pint;
        J /=pint;
        mat3d M = mat3dd(a) - Ls;
        
        double* N  = el.H (mp.m_index);
        double* Gr = el.Gr(mp.m_index);
        double* Gs = el.Gs(mp.m_index);
        
        // evaluate fluid dilatation
        double ef = GetFluidDilatation(mp, tp.alphaf);
        
        // covariant basis vectors
        vec3d gr = el.eval_deriv1(rt, mp.m_index)*m_s[iel];
        vec3d gs = el.eval_deriv2(rt, mp.m_index);
        vec3d gt = gr ^ gs;
        
        // evaluate fluid pressure
        double p = pfsi->Fluid()->Pressure(ef);

        vec3d f = gt*pfsi->Fluid()->GetElastic()->Tangent_Strain(ef, 0);

        //TODO include second order gradgrad term for surface
        
        vec3d gcnt[2], gcntp[2];
        ps->ContraBaseVectors(el, mp.m_index, gcnt);
        ps->ContraBaseVectorsP(el, mp.m_index, gcntp);
        for (int i = 0; i<neln; ++i)
            gradN[i] = ((gcnt[0] * alpha + gcntp[0] * (1 - alpha))*(Gr[i]*m_s[iel]) +
            (gcnt[1] * alpha + gcntp[1] * (1 - alpha))*Gs[i]);
        
        // calculate stiffness component
        int i = dof_a.index;
        int j = dof_b.index;
        vec3d v = gr*Gs[j] - gs*Gr[j];
        mat3d A; A.skew(v);
        mat3d Kv1 = vdotTdotv(gt, cv, gradN[j]);
        mat3d Kv2 = vdotTdotv(gt, cv, gradN[j]);
        mat3d Kw = vdotTdotv(gt, cv, (gradN[j]-gradphif*N[j]/phif)/phif);
        
        mat3d Kuu = (sv*A + Kv1*(-Dw*phis/(phif*phif)+M) + Kv2*((gradphif&gradN[j])*2.0*phis/(phif*phif*phif)+(gradN[j]&gradphif)/(phif*phif)) - Kv2*(-gradJ&gradN[j])/J*phis/(phif*phif) - ((sv*gt)&gradN[j])*phis/phif)*N[i] - A*(N[i] * p); Kuu *= -alpha;
        mat3d Kuw = Kw*N[i]; Kuw *= -alpha;
        vec3d kuJ = svJ*gt*(N[i] * N[j]) - f*(N[i] * N[j]); kuJ *= -alpha;
        
        Kab.zero();
        Kab.sub(0, 0, Kuu);
        Kab.sub(0, 3, Kuw);
        
        Kab[0][6] -= kuJ.x;
        Kab[1][6] -= kuJ.y;
        Kab[2][6] -= kuJ.z;
    });
}

//-----------------------------------------------------------------------------
void FEBiphasicFSITraction::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    
    ar & m_s;
    
    if (ar.IsShallow() == false)
    {
        if (ar.IsSaving())
        {
            int NE = (int)m_elem.size();
            ar << NE;
            for (int i = 0; i < NE; ++i)
            {
                FEElement* pe = m_elem[i];
                int nid = (pe ? pe->GetID() : -1);
                ar << nid;
            }
        }
        else
        {
            FEMesh& mesh = ar.GetFEModel().GetMesh();
            int NE, nid;
            ar >> NE;
            m_elem.resize(NE, nullptr);
            for (int i = 0; i < NE; ++i)
            {
                ar >> nid;
                if (nid != -1)
                {
                    FEElement* pe = mesh.FindElementFromID(nid);
                    assert(pe);
                    m_elem[i] = pe;
                }
            }
        }
    }
}
