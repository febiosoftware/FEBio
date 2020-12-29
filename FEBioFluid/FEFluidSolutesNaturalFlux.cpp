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

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEFluidSolutesNaturalFlux, FESurfaceLoad)
ADD_PARAMETER(m_beta, "beta");
ADD_PARAMETER(m_isol   , "solute_id");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEFluidSolutesNaturalFlux::FEFluidSolutesNaturalFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem), m_dofC(pfem)
{
    m_beta = 1.0;
    m_isol = 1;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidSolutesNaturalFlux::Init()
{
    if (m_isol == -1) return false;
    
    // set up the dof lists
    FEModel* fem = GetFEModel();
    m_dofC.Clear();
    m_dofW.Clear();

    m_dofC.AddDof(fem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), m_isol-1));
    
    m_dofW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY));

    m_dof.Clear();
    m_dof.AddDofs(m_dofW);
    m_dof.AddDofs(m_dofC);
    
    FESurface& surf = GetSurface();
    if (FESurfaceLoad::Init() == false) return false;
    
    // get the list of FluidSolutes elements connected to this interface
    int NF = surf.Elements();
    m_elem.resize(NF);
    
    for (int j = 0; j<NF; ++j)
    {
        FESurfaceElement& el = surf.Element(j);
        // extract the first of two elements on this interface
        m_elem[j] = el.m_elem[0];
        // get its material and check if FluidSolutes
        FEMaterial* pm = fem->GetMaterial(m_elem[j]->GetMatID());
        FEFluidSolutes* pfsi = dynamic_cast<FEFluidSolutes*>(pm);
        if (pfsi == nullptr) {
            pm = fem->GetMaterial(el.m_elem[1]->GetMatID());
            pfsi = dynamic_cast<FEFluidSolutes*>(pm);
            if (pfsi == nullptr) return false;
            m_elem[j] = el.m_elem[1];
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------

void FEFluidSolutesNaturalFlux::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    
    if (ar.IsShallow() == false)
    {
        ar & m_dofC & m_dofW;
    }
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

//-----------------------------------------------------------------------------
void FEFluidSolutesNaturalFlux::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    m_psurf->LoadStiffness(LS, m_dof, m_dof, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
        
        FESurfaceElement& el = *mp.SurfaceElement();
        int iel = el.m_lid;
        int neln = el.Nodes();
        
        // get the diffusivity
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEFluidSolutes* pfsm = dynamic_cast<FEFluidSolutes*>(pm);
        double d0 = pfsm->GetSolute(m_isol-1)->m_pDiff->Free_Diffusivity(mp);
        double d0p = pfsm->GetSolute(m_isol-1)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, m_isol-1);
        
        double* N  = el.H (mp.m_index);
        double* Gr = el.Gr(mp.m_index);
        double* Gs = el.Gs(mp.m_index);
        
        // tangent vectors
        vec3d rt[FEElement::MAX_NODES];
        m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);
        vec3d dxr = el.eval_deriv1(rt, mp.m_index);
        vec3d dxs = el.eval_deriv2(rt, mp.m_index);
        
        vec3d n = dxr ^ dxs;
        //double da = n.unit();
        
        double alpha = tp.alphaf;
        
        vector<vec3d> gradN(neln);
        
        // Fluid velocity
        vec3d v = FluidVelocity(mp, tp.alphaf);
        double c = EffectiveConcentration(mp, tp.alphaf);
        //vec3d v = FluidVelocityElm(mp);
        //double c = EffectiveConcentrationElm(mp);
        vec3d gradc = EffectiveCGrad(mp);
        
        vec3d gcnt[2], gcntp[2];
        m_psurf->ContraBaseVectors(el, mp.m_index, gcnt);
        m_psurf->ContraBaseVectorsP(el, mp.m_index, gcntp);
        for (int i = 0; i<neln; ++i)
            gradN[i] = (gcnt[0] * alpha + gcntp[0] * (1 - alpha))*Gr[i] +
            (gcnt[1] * alpha + gcntp[1] * (1 - alpha))*Gs[i];
        
        Kab.zero();
        
        // calculate stiffness component
        int i = dof_a.index;
        int j = dof_b.index;
        vec3d kcw = n*(N[i]*N[j]*c*tp.alphaf);
        double kcc = n*((-gradc*d0p+v)*N[j]-gradN[j]*d0)*(tp.alphaf*N[i]);
        
        Kab[3][0] -= kcw.x;
        Kab[3][1] -= kcw.y;
        Kab[3][2] -= kcw.z;
        Kab[3][3] -= kcc;
    });
}

//-----------------------------------------------------------------------------
vec3d FEFluidSolutesNaturalFlux::FluidVelocity(FESurfaceMaterialPoint& mp, double alpha)
{
    vec3d vt[FEElement::MAX_NODES];
    FESurfaceElement& el = *mp.SurfaceElement();
    int neln = el.Nodes();
    for (int j = 0; j<neln; ++j) {
        FENode& node = m_psurf->Node(el.m_lnode[j]);
        vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*alpha + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1 - alpha);
    }
    return el.eval(vt, mp.m_index);
}

//-----------------------------------------------------------------------------
double FEFluidSolutesNaturalFlux::EffectiveConcentration(FESurfaceMaterialPoint& mp, double alpha)
{
    double c = 0;
    FESurfaceElement& el = *mp.SurfaceElement();
    double* H = el.H(mp.m_index);
    int neln = el.Nodes();
    for (int j = 0; j < neln; ++j) {
        FENode& node = m_psurf->Node(el.m_lnode[j]);
        double cj = node.get(m_dofC[0])*alpha + node.get_prev(m_dofC[0])*(1.0 - alpha);
        c += cj*H[j];
    }
    return c;
}

//-----------------------------------------------------------------------------
vec3d FEFluidSolutesNaturalFlux::EffectiveCGrad(FESurfaceMaterialPoint& pt)
{
    FESurfaceElement& face = *pt.SurfaceElement();
    int iel = face.m_lid;
    
    // Get the fluid stress from the fluid-FSI element
    vec3d gradc(0.0);
    FEElement* pe = m_elem[iel];
    int nint = pe->GaussPoints();
    for (int n = 0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
        gradc += spt.m_gradc[m_isol-1];
    }
    gradc /= nint;
    return gradc;
}

//-----------------------------------------------------------------------------
double FEFluidSolutesNaturalFlux::EffectiveConcentrationElm(FESurfaceMaterialPoint& pt)
{
    FESurfaceElement& face = *pt.SurfaceElement();
    int iel = face.m_lid;
    
    // Get the fluid stress from the fluid-FSI element
    double c = 0.0;
    FEElement* pe = m_elem[iel];
    int nint = pe->GaussPoints();
    for (int n = 0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
        c += spt.m_c[m_isol-1];
    }
    c /= nint;
    return c;
}

//-----------------------------------------------------------------------------
vec3d FEFluidSolutesNaturalFlux::FluidVelocityElm(FESurfaceMaterialPoint& pt)
{
    FESurfaceElement& face = *pt.SurfaceElement();
    int iel = face.m_lid;
    
    // Get the fluid stress from the fluid-FSI element
    vec3d v(0.0);
    FEElement* pe = m_elem[iel];
    int nint = pe->GaussPoints();
    for (int n = 0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
        FEFluidMaterialPoint& fpt = *(mp.ExtractData<FEFluidMaterialPoint>());
        v += fpt.m_vft;
    }
    v /= nint;
    return v;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesNaturalFlux::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    m_psurf->LoadVector(R, m_dofC, false, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {
        
        FESurfaceElement& el = *mp.SurfaceElement();
        
        // get the diffusivity
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FEFluidSolutes* pfsm = dynamic_cast<FEFluidSolutes*>(pm);
        double d0 = pfsm->GetSolute(m_isol-1)->m_pDiff->Free_Diffusivity(mp);
        
        // tangent vectors
        vec3d rt[FEElement::MAX_NODES];
        m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);
        vec3d dxr = el.eval_deriv1(rt, mp.m_index);
        vec3d dxs = el.eval_deriv2(rt, mp.m_index);
        
        // normal and area element
        vec3d n = dxr ^ dxs;
        //double da = n.unit();
        
        // fluid velocity
        vec3d v = FluidVelocity(mp, tp.alphaf);
        double c = EffectiveConcentration(mp, tp.alphaf);
        //vec3d v = FluidVelocityElm(mp);
        //double c = EffectiveConcentrationElm(mp);
        vec3d gradc = EffectiveCGrad(mp);
        
        // force vector (change sign for inflow vs outflow)
        double f = n*(-gradc*d0+v*c)*(m_beta);
        
        double H = dof_a.shape;
        fa[0] = H * f;
    });
}
