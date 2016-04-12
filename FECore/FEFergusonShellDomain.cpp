//
//  FEFergusonShellDomain.cpp
//  FECore
//
//  Created by Gerard Ateshian on 1/29/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEFergusonShellDomain.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "FECore/DOFS.h"
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
void FEFergusonShellDomain::InitElements()
{
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FEFergusonShellElement& el = m_Elem[i];
        int n = el.GaussPoints();
        for (int j=0; j<n; ++j) el.GetMaterialPoint(j)->Init(false);
    }
}

//-----------------------------------------------------------------------------
void FEFergusonShellDomain::Reset()
{
    for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
//! calculates referential covariant basis vectors at an integration point
void FEFergusonShellDomain::CoBaseVectors0(FEFergusonShellElement& el, int n, vec3d g[3])
{
    int i;
    
    // number of nodes
    int neln = el.Nodes();
    
    assert(neln == 4);  // for now Ferguson shell elements can only be 4-node elements

    // nodal covariant basis vectors
    vec3d r0[4], gr[4], gs[4];
    vec3d t[4];
    NodalCoBaseVectors0(el, r0, gr, gs, t);
    
    // initial thickness at nodes
    double* h0 = &el.m_h0[0];
    
    // through-thickness parametric coordinate of gauss point
    double eta = el.gt(n);
    
    // get ghat and its derivatives
    vec3d ghr = vec3d(0,0,0), ghs = vec3d(0,0,0);
    vec3d ghrr = vec3d(0,0,0), ghrs = vec3d(0,0,0), ghss = vec3d(0,0,0);
    double H = 0, Hr = 0, Hs = 0;
    for (i=0; i<neln; ++i)
    {
        const double& Mi  = el.H(n)[i];
        const double& Mri = el.Hr(n)[i];
        const double& Msi = el.Hs(n)[i];
        const double& Mrri = el.Hrr(n)[i];
        const double& Mrsi = el.Hrs(n)[i];
        const double& Mssi = el.Hss(n)[i];
        const double& Pri = el.Pr(n)[i];
        const double& Psi = el.Ps(n)[i];
        const double& Prri = el.Prr(n)[i];
        const double& Prsi = el.Prs(n)[i];
        const double& Pssi = el.Pss(n)[i];
        const double& Qri = el.Qr(n)[i];
        const double& Qsi = el.Qs(n)[i];
        const double& Qrri = el.Qrr(n)[i];
        const double& Qrsi = el.Qrs(n)[i];
        const double& Qssi = el.Qss(n)[i];
        
        H  += h0[i]*Mi;     // shell initial thicknesss at integration point
        Hr += h0[i]*Mri;    // parametric derivatives of shell thickness
        Hs += h0[i]*Msi;
        
        // mid-shell covariant basis vectors at integration point
        ghr += r0[i]*Mri + gr[i]*Pri + gs[i]*Qri;
        ghs += r0[i]*Msi + gr[i]*Psi + gs[i]*Qsi;
        
        // mid-shell parametric derivatives of covariant basis vectors at integration point
        ghrr += r0[i]*Mrri + gr[i]*Prri + gs[i]*Qrri;
        ghrs += r0[i]*Mrsi + gr[i]*Prsi + gs[i]*Qrsi;
        ghss += r0[i]*Mssi + gr[i]*Pssi + gs[i]*Qssi;
    }
    
    // unit director tn at integration point
    vec3d tn = ghr ^ ghs;
    double Ln = tn.unit();
    
    // parametric derivatives of unit director at integration point
    mat3dd I(1);
    vec3d trn = (I - (tn & tn))*((ghrr ^ ghs) + (ghr ^ ghrs))/Ln;
    vec3d tsn = (I - (tn & tn))*((ghrs ^ ghs) + (ghr ^ ghss))/Ln;
    
    // evaluate covariant basis functions at current integration point
    g[0] = ghr + (tn*Hr + trn*H)*eta;
    g[1] = ghs + (tn*Hs + tsn*H)*eta;
    g[2] = tn*H/2;
}

//-----------------------------------------------------------------------------
//! calculates referential contravariant basis vectors at an integration point
void FEFergusonShellDomain::ContraBaseVectors0(FEFergusonShellElement& el, int n, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors0(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
}

//-----------------------------------------------------------------------------
double FEFergusonShellDomain::invjac0(FEFergusonShellElement& el, double Ji[3][3], int n)
{
    vec3d gcov[3];
    CoBaseVectors0(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    
    double det = J.det();
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);
    
    // calculate the inverse of the jacobian
    double deti = 1.0 / det;
    
    Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
    Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
    Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
    Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
    Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
    
    Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
    Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
    Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
    
    return det;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to reference frame
double FEFergusonShellDomain::detJ0(FEFergusonShellElement &el, int n)
{
    vec3d gcov[3];
    CoBaseVectors0(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    
    double det = J.det();
    
    return det;
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FEFergusonShellDomain::CoBaseVectors(FEFergusonShellElement& el, int n, vec3d g[3])
{
    int i;
    
    // number of nodes
    int neln = el.Nodes();
    
    assert(neln == 4);  // for now Ferguson shell elements can only be 4-node elements
    
    vec3d rt[FEElement::MAX_NODES];
    vec3d gr[4], gs[4];
    vec3d t[4];     // nodal unit directors
    double lam[4];  // nodal stretch ratios
    NodalCoBaseVectors(el, rt, gr, gs, t, lam);
    
    // initial thickness at nodes
    double* h0 = &el.m_h0[0];

    // through-thickness parametric coordinate of gauss point
    double eta = el.gt(n);

    // get ghat and its derivatives
    vec3d ghr = vec3d(0,0,0), ghs = vec3d(0,0,0);
    vec3d ghrr = vec3d(0,0,0), ghrs = vec3d(0,0,0), ghss = vec3d(0,0,0);
    double L = 0, H = 0;
    for (i=0; i<neln; ++i)
    {
        const double& Mi  = el.H(n)[i];
        const double& Mri = el.Hr(n)[i];
        const double& Msi = el.Hs(n)[i];
        const double& Mrri = el.Hrr(n)[i];
        const double& Mrsi = el.Hrs(n)[i];
        const double& Mssi = el.Hss(n)[i];
        const double& Pri = el.Pr(n)[i];
        const double& Psi = el.Ps(n)[i];
        const double& Prri = el.Prr(n)[i];
        const double& Prsi = el.Prs(n)[i];
        const double& Pssi = el.Pss(n)[i];
        const double& Qri = el.Qr(n)[i];
        const double& Qsi = el.Qs(n)[i];
        const double& Qrri = el.Qrr(n)[i];
        const double& Qrsi = el.Qrs(n)[i];
        const double& Qssi = el.Qss(n)[i];
        
        L += Mi*lam[i]; // stretch ratio at integration point
        H += Mi*h0[i];  // shell initial thicknesss at integration point
        
        // mid-shell covariant basis vectors at integration point
        ghr += rt[i]*Mri + gr[i]*Pri + gs[i]*Qri;
        ghs += rt[i]*Msi + gr[i]*Psi + gs[i]*Qsi;

        // mid-shell parametric derivatives of covariant basis vectors at integration point
        ghrr += rt[i]*Mrri + gr[i]*Prri + gs[i]*Qrri;
        ghrs += rt[i]*Mrsi + gr[i]*Prsi + gs[i]*Qrsi;
        ghss += rt[i]*Mssi + gr[i]*Pssi + gs[i]*Qssi;
    }
    
    // unit director tn at integration point
    vec3d tn = ghr ^ ghs;
    double Ln = tn.unit();
    
    // parametric derivatives of unit director at integration point
    mat3dd I(1);
    vec3d trn = (I - (tn & tn))*((ghrr ^ ghs) + (ghr ^ ghrs))/Ln;
    vec3d tsn = (I - (tn & tn))*((ghrs ^ ghs) + (ghr ^ ghss))/Ln;
    
    // current thickness h at integration point
    double h = L*H;
    
    // parametric derivatives of h at integration point
    double hr = 0, hs = 0;
    for (i=0; i<neln; ++i)
    {
        const double& Mri = el.Hr(n)[i];
        const double& Msi = el.Hs(n)[i];
        hr += (L*h0[i] + lam[i]*H)*Mri;
        hs += (L*h0[i] + lam[i]*H)*Msi;
    }
    
    // evaluate covariant basis functions at current integration point
    g[0] = ghr + (tn*hr + trn*h)*eta;
    g[1] = ghs + (tn*hs + tsn*h)*eta;
    g[2] = tn*h/2;
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
//! and returns results of intermediate calculations
void FEFergusonShellDomain::MidShellSurface(FEFergusonShellElement& el, int n,
                                            double& L, double& Lr, double& Ls,
                                            double& H, double& Hr, double& Hs,
                                            vec3d& ghr, vec3d& ghs,
                                            vec3d& ghrr, vec3d& ghrs, vec3d& ghss,
                                            vec3d& tn,vec3d& trn, vec3d& tsn)
{
    int i;
    
    // number of nodes
    int neln = el.Nodes();
    
    assert(neln == 4);  // for now Ferguson shell elements can only be 4-node elements
    
    vec3d rt[FEElement::MAX_NODES];
    vec3d gr[4], gs[4];
    vec3d t[4];     // nodal unit directors
    double lam[4];  // nodal stretch ratios
    NodalCoBaseVectors(el, rt, gr, gs, t, lam);
    
    // initial thickness at nodes
    double* h0 = &el.m_h0[0];
    
    // get ghat and its derivatives
    ghr = vec3d(0,0,0), ghs = vec3d(0,0,0);
    ghrr = vec3d(0,0,0), ghrs = vec3d(0,0,0), ghss = vec3d(0,0,0);
    L = 0; Lr = 0; Ls = 0;
    H = 0; Hr = 0; Hs = 0;
    for (i=0; i<neln; ++i)
    {
        const double& Mi  = el.H(n)[i];
        const double& Mri = el.Hr(n)[i];
        const double& Msi = el.Hs(n)[i];
        const double& Mrri = el.Hrr(n)[i];
        const double& Mrsi = el.Hrs(n)[i];
        const double& Mssi = el.Hss(n)[i];
        const double& Pri = el.Pr(n)[i];
        const double& Psi = el.Ps(n)[i];
        const double& Prri = el.Prr(n)[i];
        const double& Prsi = el.Prs(n)[i];
        const double& Pssi = el.Pss(n)[i];
        const double& Qri = el.Qr(n)[i];
        const double& Qsi = el.Qs(n)[i];
        const double& Qrri = el.Qrr(n)[i];
        const double& Qrsi = el.Qrs(n)[i];
        const double& Qssi = el.Qss(n)[i];
        
        L += Mi*lam[i]; // stretch ratio at integration point
        Lr += Mri*lam[i];
        Ls += Msi*lam[i];
        H += Mi*h0[i];  // shell initial thicknesss at integration point
        Hr += Mri*h0[i];
        Hs += Msi*h0[i];
        
        // mid-shell covariant basis vectors at integration point
        ghr += rt[i]*Mri + gr[i]*Pri + gs[i]*Qri;
        ghs += rt[i]*Msi + gr[i]*Psi + gs[i]*Qsi;
        
        // mid-shell parametric derivatives of covariant basis vectors at integration point
        ghrr += rt[i]*Mrri + gr[i]*Prri + gs[i]*Qrri;
        ghrs += rt[i]*Mrsi + gr[i]*Prsi + gs[i]*Qrsi;
        ghss += rt[i]*Mssi + gr[i]*Pssi + gs[i]*Qssi;
    }
    
    // unit director tn at integration point
    tn = ghr ^ ghs;
    double Ln = tn.unit();
    
    // parametric derivatives of unit director at integration point
    mat3dd I(1);
    trn = (I - (tn & tn))*((ghrr ^ ghs) + (ghr ^ ghrs))/Ln;
    tsn = (I - (tn & tn))*((ghrs ^ ghs) + (ghr ^ ghss))/Ln;
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FEFergusonShellDomain::ContraBaseVectors(FEFergusonShellElement& el, int n, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the current frame at
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
double FEFergusonShellDomain::invjact(FEFergusonShellElement& el, double Ji[3][3], int n)
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    
    double det = J.det();
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);
    
    // calculate the inverse of the jacobian
    double deti = 1.0 / det;
    
    Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
    Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
    Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
    Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
    Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
    
    Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
    Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
    Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
    
    return det;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to reference frame
double FEFergusonShellDomain::detJt(FEFergusonShellElement &el, int n)
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    
    double det = J.det();
    
    return det;
}

//-----------------------------------------------------------------------------
double FEFergusonShellDomain::defgrad(FEFergusonShellElement& el, mat3d& F, int n)
{
    vec3d gcov[3], Gcnt[3];
    CoBaseVectors(el, n, gcov);
    ContraBaseVectors0(el, n, Gcnt);

    F = (gcov[0] & Gcnt[0]) + (gcov[1] & Gcnt[1]) + (gcov[2] & Gcnt[2]);
    
    double V = F.det();
    if (V <= 0) throw NegativeJacobian(el.GetID(), n, V, &el);
    
    return V;
}

//-----------------------------------------------------------------------------
//! calculates nodal covariant basis vectors
void FEFergusonShellDomain::NodalCoBaseVectors0(FEFergusonShellElement& el, vec3d* r0, vec3d* gr, vec3d* gs, vec3d* t)
{
    int i;
    
    // number of nodes
    int neln = el.Nodes();
    
    assert(neln == 4);  // for now Ferguson shell elements can only be 4-node elements
    
    // initial nodal coordinates and directors
    vec3d D0[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r0[i] = ni.m_r0;
        D0[i] = el.m_D0[i];
        t[i] = D0[i];
        t[i].unit();
    }
    
    // evaluate four corner tangents assuming bilinear element
    vec3d gtr[4], gts[4];
    gtr[0] = gtr[1] = (r0[1] - r0[0])/2;
    gtr[2] = gtr[3] = (r0[2] - r0[3])/2;
    gts[0] = gts[3] = (r0[3] - r0[0])/2;
    gts[1] = gts[2] = (r0[2] - r0[1])/2;
    
    // project each corner tangent on plane normal to director t
    mat3dd I(1);
    for (i=0; i<4; ++i) {
        gr[i] = (I - (t[i] & t[i]))*gtr[i];
        gs[i] = (I - (t[i] & t[i]))*gts[i];
    }
}

//-----------------------------------------------------------------------------
//! calculates nodal covariant basis vectors
void FEFergusonShellDomain::NodalCoBaseVectors(FEFergusonShellElement& el, vec3d* rt, vec3d* gr, vec3d* gs, vec3d* t, double* lam)
{
    int i;
    
    // number of nodes
    int neln = el.Nodes();
    
    assert(neln == 4);  // for now Ferguson shell elements can only be 4-node elements
    
    // initial nodal coordinates and directors
    vec3d Dt[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        rt[i] = ni.m_rt;
        Dt[i] = el.m_D0[i] + ni.get_vec3d(m_dofU, m_dofV, m_dofW);
        t[i] = Dt[i];
        lam[i] = t[i].unit();
    }
    
    // evaluate four corner tangents assuming bilinear element
    vec3d gtr[4], gts[4];
    gtr[0] = gtr[1] = (rt[1] - rt[0])/2;
    gtr[2] = gtr[3] = (rt[2] - rt[3])/2;
    gts[0] = gts[3] = (rt[3] - rt[0])/2;
    gts[1] = gts[2] = (rt[2] - rt[1])/2;
    
    // project each corner tangent on plane normal to director t
    mat3dd I(1);
    for (i=0; i<neln; ++i) {
        mat3ds Imtt = I - dyad(t[i]);
        gr[i] = Imtt*gtr[i];
        gs[i] = Imtt*gts[i];
    }
}

//-----------------------------------------------------------------------------
//! calculates nodal covariant basis torsional tangents in current frame
void FEFergusonShellDomain::NodalCoBaseTorsionTangents(FEFergusonShellElement& el, vec3d* rt, vec3d* t, double* lam, mat3d* Gr, mat3d* Gs)
{
    int i;
    
    // number of nodes
    int neln = el.Nodes();
    
    assert(neln == 4);  // for now Ferguson shell elements can only be 4-node elements
    
    // evaluate four corner tangents assuming bilinear element
    vec3d gtr[4], gts[4];
    gtr[0] = gtr[1] = (rt[1] - rt[0])/2;
    gtr[2] = gtr[3] = (rt[2] - rt[3])/2;
    gts[0] = gts[3] = (rt[3] - rt[0])/2;
    gts[1] = gts[2] = (rt[2] - rt[1])/2;
    
    for (i=0; i<neln; ++i) {
        mat3ds Imtt = mat3dd(1) - dyad(t[i]);
        Gr[i] = (mat3dd(t[i]*gtr[i]) + (t[i] & gtr[i]))*Imtt/(-lam[i]);
        Gs[i] = (mat3dd(t[i]*gts[i]) + (t[i] & gts[i]))*Imtt/(-lam[i]);
    }
}

//-----------------------------------------------------------------------------
void FEFergusonShellDomain::Serialize(DumpStream &ar)
{
    if (ar.IsShallow())
    {
        int NEL = (int) m_Elem.size();
        for (int i=0; i<NEL; ++i)
        {
            FEFergusonShellElement& el = m_Elem[i];
            int nint = el.GaussPoints();
            for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
        }
    }
    else
    {
        if (ar.IsSaving())
        {
            ar << m_Node;
            
            for (size_t i=0; i<m_Elem.size(); ++i)
            {
                FEFergusonShellElement& el = m_Elem[i];
                ar << el.Type();
                
                ar << el.GetMatID();
                ar << el.GetID();
                ar << el.m_node;
                
                ar << el.m_h0;
                ar << el.m_D0;
                
                for (int j=0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
            }
        }
        else
        {
            int n, mat, nid;
            
            ar >> m_Node;
            
            FEModel& fem = ar.GetFEModel();
            FEMaterial* pmat = GetMaterial();
            assert(pmat);
            
            for (size_t i=0; i<m_Elem.size(); ++i)
            {
                FEFergusonShellElement& el = m_Elem[i];
                ar >> n;
                
                el.SetType(n);
                
                ar >> mat; el.SetMatID(mat);
                ar >> nid; el.SetID(nid);
                ar >> el.m_node;
                
                ar >> el.m_h0;
                ar >> el.m_D0;
                
                for (int j=0; j<el.GaussPoints(); ++j)
                {
                    el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
                    el.GetMaterialPoint(j)->Serialize(ar);
                }
            }
        }
    }
}
