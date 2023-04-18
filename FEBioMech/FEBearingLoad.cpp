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
#include "FEBearingLoad.h"
#include "FEBioMech.h"
#include <FECore/FESurface.h>
#include <FECore/FEFacetSet.h>
#include <FECore/FEMesh.h>
#include <FECore/QuadricFit.h>
#include <FECore/log.h>
#include <string>

//-----------------------------------------------------------------------------
// Parameter block for bearing loads
BEGIN_FECORE_CLASS(FEBearingLoad, FESurfaceLoad)
    ADD_PARAMETER(m_scale   , "scale");
    ADD_PARAMETER(m_force   , "force");
    ADD_PARAMETER(m_bsymm   , "symmetric_stiffness");
    ADD_PARAMETER(m_blinear , "linear");
    ADD_PARAMETER(m_bshellb , "shell_bottom");
    ADD_PARAMETER(m_profile , "profile");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEBearingLoad::FEBearingLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_scale = 1.0;
    m_force = vec3d(0,0,0);
    m_bsymm = true;
    m_bshellb = false;
    m_blinear = false;
    m_profile = P_SINE;
}

//-----------------------------------------------------------------------------
bool FEBearingLoad::Init()
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
    
    std::vector<vec3d> pc;
    
    FESurface& surf = GetSurface();
    surf.SetShellBottom(m_bshellb);
    FENodeList nl = surf.GetNodeList();
    FEMesh* mesh = surf.GetMesh();
    
    // figure out the pressure load distribution on the bearing surface.
    // 1. Fit a quadric to the bearing surface
    // 2. Confirm that the quadric is a cylinder
    // 3. Set the cylinder axis ez and create local coordinate system
    //    based on orientation of bearing force (er)
    // 4. For each face of bearing surface, get cosine of angle theta of its normal relative to er (cq)
    // 5. Use cq to evaluate the pressure profile and projected area, scale to match prescribed force
    
    // collect all the surface nodal positions into pc
    if (m_bshellb == false) {
        for (int i=0; i<nl.Size(); ++i)
            pc.push_back(mesh->Node(nl[i]).m_r0);
    }
    else {
        for (int i=0; i<nl.Size(); ++i)
            pc.push_back(mesh->Node(nl[i]).s0());
    }
    
    // find the best fit quadric
    QuadricFit* fit = new QuadricFit();
    fit->Fit(pc);
    
    QuadricFit::Q_TYPE qtype = fit->GetType();

    bool cylindrical = qtype & (QuadricFit::Q_CIRCULAR_CYLINDER |
                                QuadricFit::Q_ELLIPTIC_CYLINDER |
                                QuadricFit::Q_PARABOLIC_CYLINDER|
                                QuadricFit::Q_HYPERBOLIC_CYLINDER);
    
    if (cylindrical == false) {
        std::string quadric = fit->GetStringType(qtype);
        feLogError("Bearing load surface is %s, not cylindrical!\n",quadric.c_str());
        return false;
    }
    
    // extract cylinder parameters
//    vec3d c = fit->m_rc;                // cylinder origin
//    double R = sqrt(1./fit->m_c2.x);    // radius
    vec3d ez = fit->m_ax[2];            // cylinder axis
    vec3d et = (ez ^ m_force).normalized();
    m_er = (et ^ ez).normalized();  // radial force direction
    
    // evaluate dot product of bearing force and axis
    double dp = m_force*ez;
    double eps = 1e-3;
    if (fabs(dp) > m_force.norm()*eps) {
        feLogWarning("Axial component of bearing force (%g) is ignored!\n",dp);
    }
    
    // create surface map for prescribed pressure
    m_pc = new FESurfaceMap(FE_DOUBLE);
    m_pc->Create(surf.GetFacetSet(), 0.0, FMT_MULT);
    int m = m_pc->MaxNodes();

    // evaluate projected area, which is needed to evaluate peak pressure
    double A = 0;       // effective projected area
    for (int i=0; i<surf.Elements(); ++i) {
        // get the surface element
        FESurfaceElement& el = surf.Element(i);
        // get (and negate inward) normal at face centroid
        vec3d n = -surf.SurfaceNormal(el, el.cr(), el.cs());
        // get face area
        double area = surf.FaceArea(el);
        // cosine of angle between force and surface normal
        double cq = n*m_er;
        // if cq is positive
        if (cq > 0) {
            // face projected area
            double parea = area*cq;
            // add to effective projected area
            switch (m_profile) {
                case P_SINE: A += cq*parea; break;
                case P_PARA: A += pow(cq,2)*parea; break;
                default: break;
            }
        }
    }
    // evaluate peak pressure p0 for a unit force magnitude
    double p0 = 1.0/A;
    
    // prescribe bearing pressure on each face
    for (int i=0; i<surf.Elements(); ++i) {
        // get the surface element
        FESurfaceElement& el = surf.Element(i);
        // get (and negate inward) normal at face centroid
        vec3d n = -surf.SurfaceNormal(el, el.cr(), el.cs());
        // cosine of angle between force and surface normal
        double cq = n*m_er;
        // if cq is positive
        if (cq > 0) {
            double p = p0;
            switch (m_profile) {
                case P_SINE: p *= cq; break;
                case P_PARA: p *= pow(cq,2); break;
                default: p = 0; break;
            }
            for (int j = 0; j < m; ++j)
                m_pc->setValue(el.m_lid, j, p);
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! update projected area of bearing surface when nonlinear flag is set
void FEBearingLoad::Update()
{
    if (m_blinear == false) {
        FESurface& surf = GetSurface();
        surf.SetShellBottom(m_bshellb);
        int m = m_pc->MaxNodes();

        // evaluate projected area, which is needed to evaluate peak pressure
        double A = 0;       // effective projected area
        for (int i=0; i<surf.Elements(); ++i) {
            // get the surface element
            FESurfaceElement& el = surf.Element(i);
            // get (and negate inward) normal at face centroid
            vec3d n = -surf.SurfaceNormal(el, el.cr(), el.cs());
            // get face area
            double area = surf.FaceArea(el);
            // cosine of angle between force and surface normal
            double cq = n*m_er;
            // if cq is positive
            if (cq > 0) {
                // face projected area
                double parea = area*cq;
                // add to effective projected area
                switch (m_profile) {
                    case P_SINE: A += cq*parea; break;
                    case P_PARA: A += pow(cq,2)*parea; break;
                    default: break;
                }
            }
        }
        // evaluate peak pressure p0 for a unit force magnitude
        double p0 = 1.0/A;
        
        // prescribe bearing pressure on each face
        for (int i=0; i<surf.Elements(); ++i) {
            // get the surface element
            FESurfaceElement& el = surf.Element(i);
            // get (and negate inward) normal at face centroid
            vec3d n = -surf.SurfaceNormal(el, el.cr(), el.cs());
            // cosine of angle between force and surface normal
            double cq = n*m_er;
            // if cq is positive
            if (cq > 0) {
                double p = p0;
                switch (m_profile) {
                    case P_SINE: p *= cq; break;
                    case P_PARA: p *= pow(cq,2); break;
                    default: p = 0; break;
                }
                for (int j = 0; j < m; ++j)
                    m_pc->setValue(el.m_lid, j, p);
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEBearingLoad::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow() == false)
    {
        ar & m_er;

        if (ar.IsSaving())
            m_pc->Serialize(ar);
        else
        {
            m_pc = new FESurfaceMap(FE_DOUBLE);
            m_pc->Serialize(ar);
        }
    }
}

//-----------------------------------------------------------------------------
//! evaluate bearing pressure
double FEBearingLoad::ScalarLoad(FESurfaceMaterialPoint& mp)
{
    // evaluate pressure at this material point
    double P = m_pc->value(mp)*m_scale*m_force.norm();
    return P;
}

//-----------------------------------------------------------------------------
void FEBearingLoad::LoadVector(FEGlobalVector& R)
{
    FESurface& surf = GetSurface();
    surf.SetShellBottom(m_bshellb);
    
    // evaluate the integral
    surf.LoadVector(R, m_dof, m_blinear, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {
        
        // evaluate pressure at this material point
        double P = -ScalarLoad(pt);
        if (m_bshellb) P = -P;

        double J = (pt.dxr ^ pt.dxs).norm();
        
        // force vector
        vec3d N = (pt.dxr ^ pt.dxs); N.unit();
        vec3d t = N*P;
        
        double H_u = dof_a.shape;
        
        val[0] = H_u*t.x*J;
        val[1] = H_u*t.y*J;
        val[2] = H_u*t.z*J;
    });
}

//-----------------------------------------------------------------------------
void FEBearingLoad::StiffnessMatrix(FELinearSystem& LS)
{
    // Don't calculate stiffness for a linear load
    if (m_blinear) return;
    
    FESurface& surf = GetSurface();
    surf.SetShellBottom(m_bshellb);
    
    // evaluate the integral
    surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& kab) {
        
        // evaluate pressure at this material point
        double P = -ScalarLoad(mp);
        if (m_bshellb) P = -P;

        double H_i  = dof_a.shape;
        double Gr_i = dof_a.shape_deriv_r;
        double Gs_i = dof_a.shape_deriv_s;
        
        double H_j  = dof_b.shape;
        double Gr_j = dof_b.shape_deriv_r;
        double Gs_j = dof_b.shape_deriv_s;
        
        vec3d vab(0,0,0);
        if (m_bsymm)
            vab = (mp.dxr*(H_j * Gs_i - H_i * Gs_j) - mp.dxs*(H_j * Gr_i - H_i * Gr_j)) * 0.5*P;
        else
            vab = (mp.dxs*Gr_j - mp.dxr*Gs_j)*(P*H_i);
        
        mat3da K(vab);
        kab.set(0, 0, K);
    });
}
