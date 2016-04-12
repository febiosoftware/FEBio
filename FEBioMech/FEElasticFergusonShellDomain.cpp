//
//  FEElasticFergusonShellDomain.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 1/29/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEElasticFergusonShellDomain.h"
#include "FEElasticMaterial.h"
#include "FEBodyForce.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <math.h>

//-----------------------------------------------------------------------------
FEElasticFergusonShellDomain::FEElasticFergusonShellDomain(FEModel* pfem) : FEFergusonShellDomain(&pfem->GetMesh()), FEElasticDomain(pfem)
{
    m_pMat = 0;
    m_dofU = pfem->GetDOFIndex("u");
    m_dofV = pfem->GetDOFIndex("v");
    m_dofW = pfem->GetDOFIndex("w");
}

//-----------------------------------------------------------------------------
void FEElasticFergusonShellDomain::SetMaterial(FEMaterial* pmat)
{
    m_pMat = dynamic_cast<FESolidMaterial*>(pmat);
}

//-----------------------------------------------------------------------------
bool FEElasticFergusonShellDomain::Initialize(FEModel& mdl)
{
    // initialize base class
    FEFergusonShellDomain::Initialize(mdl);
    
    // error flag (set true on error)
    bool bmerr = false;
    
    // set the local coordinate system for each integration point
    FECoordSysMap* pmap = m_pMat->GetElasticMaterial()->GetCoordinateSystemMap();
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        // unpack element data
        FEFergusonShellElement& el = m_Elem[i];
        
        // set the local element coordinates
        if (pmap)
        {
            for (int n=0; n<el.GaussPoints(); ++n)
            {
                FEElasticMaterialPoint& pt = *el.GetMaterialPoint(n)->ExtractData<FEElasticMaterialPoint>();
                pt.m_Q = pmap->LocalElementCoord(el, n);
            }
        }
    }
    
    // check for initially inverted shells
    for (int i=0; i<Elements(); ++i)
    {
        FEFergusonShellElement& el = Element(i);
        int nint = el.GaussPoints();
        for (int n=0; n<nint; ++n)
        {
            double J0 = detJ0(el, n);
            if (J0 <= 0)
            {
                felog.printf("**************************** E R R O R ****************************\n");
                felog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.GetID());
                felog.printf("Jacobian = %lg\n", J0);
                felog.printf("Did you use the right node numbering?\n");
                felog.printf("Nodes:");
                for (int l=0; l<el.Nodes(); ++l)
                {
                    felog.printf("%d", el.m_node[l]+1);
                    if (l+1 != el.Nodes()) felog.printf(","); else felog.printf("\n");
                }
                felog.printf("*******************************************************************\n\n");
                bmerr = true;
            }
        }
    }
    
    return (bmerr == false);
}

//-----------------------------------------------------------------------------
void FEElasticFergusonShellDomain::Activate()
{
    for (int i=0; i<Nodes(); ++i)
    {
        FENode& node = Node(i);
        if (node.m_bexclude == false)
        {
            if (node.m_rid < 0)
            {
                node.m_ID[m_dofX] = DOF_ACTIVE;
                node.m_ID[m_dofY] = DOF_ACTIVE;
                node.m_ID[m_dofZ] = DOF_ACTIVE;
            }
            
            if (node.m_bshell)
            {
                node.m_ID[m_dofU] = DOF_ACTIVE;
                node.m_ID[m_dofV] = DOF_ACTIVE;
                node.m_ID[m_dofW] = DOF_ACTIVE;
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Calculates the forces due to the stress
void FEElasticFergusonShellDomain::InternalForces(FEGlobalVector& R)
{
    // element force vector
    vector<double> fe;
    
    vector<int> lm;
    
    int NS = m_Elem.size();
    for (int i=0; i<NS; ++i)
    {
        // get the element
        FEFergusonShellElement& el = m_Elem[i];
        
        // create the element force vector and initialize to zero
        int ndof = 6*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate element's internal force
        ElementInternalForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble the residual
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for shell elements
//! Note that we use a one-point gauss integration rule for the thickness
//! integration. This will integrate linear functions exactly.

void FEElasticFergusonShellDomain::ElementInternalForce(FEFergusonShellElement& el, vector<double>& fe)
{
    int i, j, n;
    double fr[4][4] = {
        {-0.5, 0.5, 0.0, 0.0},
        {-0.5, 0.5, 0.0, 0.0},
        { 0.0, 0.0, 0.5,-0.5},
        { 0.0, 0.0, 0.4,-0.5}
    };
    double fs[4][4] = {
        {-0.5, 0.0, 0.0, 0.5},
        { 0.0,-0.5, 0.5, 0.0},
        { 0.0,-0.5, 0.5, 0.0},
        {-0.5, 0.0, 0.0, 0.5}
    };
    
    // nodal covariant basis vectors and associated data
    double lam[4];
    vec3d rt[4], gr[4], gs[4], t[4];
    mat3d Gr[4], Gs[4];
    NodalCoBaseVectors(el, rt, gr, gs, t, lam);
    NodalCoBaseTorsionTangents(el, rt, t, lam, Gr, Gs);
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double detJ;
    
    const double *M, *Mr, *Ms;
    const double *Mrr, *Mrs, *Mss;
    const double *P, *Pr, *Ps;
    const double *Prr, *Prs, *Pss;
    const double *Q, *Qr, *Qs;
    const double *Qrr, *Qrs, *Qss;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*	gw = el.GaussWeights();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEElasticMaterialPoint& pt = *(el.GetMaterialPoint(n)->ExtractData<FEElasticMaterialPoint>());
        
        // calculate the jacobian
        detJ = detJt(el, n);
        
        detJ *= gw[n];
        
        // get the stress vector for this integration point
        mat3ds& T = pt.m_s;
        
        // get relevant data for this integration point
        vec3d gcnt[3];
        vec3d ghr, ghs, ghrr, ghrs, ghss;
        vec3d tn, trn, tsn;
        double H, Hr, Hs, L, Lr, Ls;
        ContraBaseVectors(el, n, gcnt);
        MidShellSurface(el, n, L, Lr, Ls, H, Hr, Hs, ghr, ghs, ghrr, ghrs, ghss, tn, trn, tsn);
        double eta = el.gt(n);
        double gmag = (ghr ^ ghs).norm();
        mat3ds Imttg = (mat3dd(1) - dyad(tn))/gmag;
        mat3ds Imttgr = (dyads(tn, trn) + Imttg*(tn*((ghrr ^ ghs) + (ghr ^ ghrs))))/(-gmag);
        mat3ds Imttgs = (dyads(tn, tsn) + Imttg*(tn*((ghrs ^ ghs) + (ghr ^ ghss))))/(-gmag);
        mat3d Aghr, Aghs, Aghrr, Aghrs, Aghss;
        Aghr.skew(ghr); Aghs.skew(ghs);
        Aghrr.skew(ghrr); Aghrs.skew(ghrs); Aghss.skew(ghss);
        
        M = el.H(n); Mr = el.Hr(n); Ms = el.Hs(n);
        Mrr = el.Hrr(n); Mrs = el.Hrs(n); Mss = el.Hss(n);
        P = el.P(n); Pr = el.Pr(n); Ps = el.Ps(n);
        Prr = el.Prr(n); Prs = el.Prs(n); Pss = el.Pss(n);
        Q = el.Q(n); Qr = el.Qr(n); Qs = el.Qs(n);
        Qrr = el.Qrr(n); Qrs = el.Qrs(n); Qss = el.Qss(n);
        
        mat3ds Khr[4], Khs[4];
        mat3ds Khrr[4], Khrs[4], Khss[4];
        mat3d Qhr[4], Qhs[4];
        mat3d Qhrr[4], Qhrs[4], Qhss[4];
        mat3d Kr[4], Ks[4], Kt[4];
        mat3d dQr[4], dQs[4], dQt[4];
        
        mat3dd I(1);
        
        for (i=0; i<neln; ++i)
        {
            mat3ds pr, ps, prr, prs, pss;
            pr.zero(); ps.zero(); prr.zero(); prs.zero(); pss.zero();
            for (j=0; j<neln; ++j)
            {
                mat3ds Imtt = I - dyad(t[j]);
                pr += Imtt*(Pr[j]*fr[j][i] + Qr[j]*fs[j][i]);
                ps += Imtt*(Ps[j]*fr[j][i] + Qs[j]*fs[j][i]);
                prr += Imtt*(Prr[j]*fr[j][i] + Qrr[j]*fs[j][i]);
                prs += Imtt*(Prs[j]*fr[j][i] + Qrs[j]*fs[j][i]);
                pss += Imtt*(Pss[j]*fr[j][i] + Qss[j]*fs[j][i]);
            }
            
            Khr[i] = mat3dd(Mr[i]) + pr;
            Khs[i] = mat3dd(Ms[i]) + ps;
            Khrr[i] = mat3dd(Mrr[i]) + prr;
            Khrs[i] = mat3dd(Mrs[i]) + prs;
            Khss[i] = mat3dd(Mss[i]) + pss;
            Qhr[i] = Gr[i]*Pr[i] + Gs[i]*Qr[i];
            Qhs[i] = Gr[i]*Ps[i] + Gs[i]*Qs[i];
            Qhrr[i] = Gr[i]*Prr[i] + Gs[i]*Qrr[i];
            Qhrs[i] = Gr[i]*Prs[i] + Gs[i]*Qrs[i];
            Qhss[i] = Gr[i]*Pss[i] + Gs[i]*Qss[i];
            
            mat3d AK = Imttg*(Aghr*Khs[i] - Aghs*Khr[i]);
            Kr[i] = Khr[i] + AK*((L*Hr + Lr*H)*eta/2)
            + Imttgr*(Aghr*Khs[i] - Aghs*Khr[i])*(L*H*eta/2)
            + Imttg*(Aghrr*Khs[i] - Aghrs*Khr[i] + Aghr*Khrs[i] - Aghs*Khrr[i])*(L*H*eta/2);
            Ks[i] = Khs[i] + AK*((L*Hs + Ls*H)*eta/2)
            + Imttgs*(Aghr*Khs[i] - Aghs*Khr[i])*(L*H*eta/2)
            + Imttg*(Aghrs*Khs[i] - Aghss*Khr[i] + Aghr*Khss[i] - Aghs*Khrs[i])*(L*H*eta/2);
            Kt[i] = AK*(L*H/2);
            
            mat3d AQ = (tn & t[i])*M[i] + Imttg*(Aghr*Qhs[i] - Aghs*Qhr[i])*L;
            dQr[i] = Qhr[i] + AQ*(Hr*eta/2)
            + (((tn*Mr[i] + trn*M[i]) & t[i]) + (Imttg*Lr + Imttgr*L)*(Aghr*Qhs[i] - Aghs*Qhr[i]))*(H*eta/2)
            + Imttg*(Aghrr*Qhs[i] - Aghrs*Qhr[i] + Aghr*Qhrs[i] - Aghs*Qhrr[i])*(L*H*eta/2);
            dQs[i] = Qhs[i] + AQ*(Hs*eta/2)
            + (((tn*Ms[i] + tsn*M[i]) & t[i]) + (Imttg*Ls + Imttgs*L)*(Aghr*Qhs[i] - Aghs*Qhr[i]))*(H*eta/2)
            + Imttg*(Aghrs*Qhs[i] - Aghss*Qhr[i] + Aghr*Qhss[i] - Aghs*Qhrs[i])*(L*H*eta/2);
            dQt[i] = AQ*(H/2);
            
            vec3d fu =  Kr[i].transpose()*T*gcnt[0] +  Ks[i].transpose()*T*gcnt[1] +  Kt[i].transpose()*T*gcnt[2];
            vec3d fd = dQr[i].transpose()*T*gcnt[0] + dQs[i].transpose()*T*gcnt[1] + dQt[i].transpose()*T*gcnt[2];
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector

            fe[6*i  ] -= fu.x*detJ;
            fe[6*i+1] -= fu.y*detJ;
            fe[6*i+2] -= fu.z*detJ;
            
            fe[6*i+3] -= fd.x*detJ;
            fe[6*i+4] -= fd.y*detJ;
            fe[6*i+5] -= fd.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
void FEElasticFergusonShellDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    // element force vector
    vector<double> fe;
    
    vector<int> lm;
    
    int NS = m_Elem.size();
    for (int i=0; i<NS; ++i)
    {
        // get the element
        FEFergusonShellElement& el = m_Elem[i];
        
        // create the element force vector and initialize to zero
        int ndof = 6*el.Nodes();
        fe.assign(ndof, 0);
        
        // apply body forces to shells
        ElementBodyForce(BF, el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble the residual
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! Calculates element body forces for shells

void FEElasticFergusonShellDomain::ElementBodyForce(FEBodyForce& BF, FEFergusonShellElement& el, vector<double>& fe)
{
    double dens = m_pMat->Density();
    
    // calculate the average thickness
    double* h0 = &el.m_h0[0], gt, za;
    
    // integration weights
    double* gw = el.GaussWeights();
    double *Hn, detJ;
    
    // loop over integration points
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    // nodal coordinates
    vec3d r0[4], rt[4];
    for (int i=0; i<neln; ++i)
    {
        r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
        rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
    }
    
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
        pt.m_r0 = el.Evaluate(r0, n);
        pt.m_rt = el.Evaluate(rt, n);
        
        detJ = detJ0(el, n)*gw[n];
        Hn  = el.H(n);
        gt = el.gt(n);
        
        // get the force
        vec3d f = BF.force(mp);
        
        for (int i=0; i<neln; ++i)
        {
            za = 0.5*gt*h0[i];
            
            fe[6*i  ] -= Hn[i]*f.x*dens*detJ;
            fe[6*i+1] -= Hn[i]*f.y*dens*detJ;
            fe[6*i+2] -= Hn[i]*f.z*dens*detJ;
            
            fe[6*i+3] -= za*Hn[i]*dens*f.x*detJ;
            fe[6*i+4] -= za*Hn[i]*dens*f.y*detJ;
            fe[6*i+5] -= za*Hn[i]*dens*f.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------

void FEElasticFergusonShellDomain::StiffnessMatrix(FESolver* psolver)
{
    FEModel& fem = psolver->GetFEModel();
    
    matrix ke;
    
    vector<int> lm;
    
    int NS = m_Elem.size();
    for (int iel=0; iel<NS; ++iel)
    {
        FEFergusonShellElement& el = m_Elem[iel];
        
        // get the elements material
        FEMaterial* pmat = m_pMat;
        
        // create the element's stiffness matrix
        int ndof = 6*el.Nodes();
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementStiffness(iel, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
        
        if (fem.GetCurrentStep()->GetPrintLevel() == FE_PRINT_MINOR_ITRS_EXP)
        {
            fprintf(stderr, "Calculating stiffness matrix: %.1lf %% \r", 100.0*(iel)/ NS);
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates the shell element stiffness matrix

void FEElasticFergusonShellDomain::ElementStiffness(int iel, matrix& ke)
{
    FEFergusonShellElement& el = Element(iel);
    
    int i, i6, j, j6, n;
    
    double fr[4][4] = {
        {-0.5, 0.5, 0.0, 0.0},
        {-0.5, 0.5, 0.0, 0.0},
        { 0.0, 0.0, 0.5,-0.5},
        { 0.0, 0.0, 0.4,-0.5}
    };
    double fs[4][4] = {
        {-0.5, 0.0, 0.0, 0.5},
        { 0.0,-0.5, 0.5, 0.0},
        { 0.0,-0.5, 0.5, 0.0},
        {-0.5, 0.0, 0.0, 0.5}
    };
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    const int ndof = 6*neln;
    
    // nodal covariant basis vectors and associated data
    double lam[4];
    vec3d rt[4], gr[4], gs[4], t[4];
    mat3d Gr[4], Gs[4];
    NodalCoBaseVectors(el, rt, gr, gs, t, lam);
    NodalCoBaseTorsionTangents(el, rt, t, lam, Gr, Gs);
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double detJ;
    
    const double *M, *Mr, *Ms;
    const double *Mrr, *Mrs, *Mss;
    const double *P, *Pr, *Ps;
    const double *Prr, *Prs, *Pss;
    const double *Q, *Qr, *Qs;
    const double *Qrr, *Qrs, *Qss;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    // calculate element stiffness matrix
    ke.zero();
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *(el.GetMaterialPoint(n));
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        
        // calculate the jacobian
        detJ = detJt(el, n);
        
        detJ *= gw[n];
        
        // get the stress and elasticity tensors for this integration point
        mat3ds T = pt.m_s;
        tens4ds C = m_pMat->Tangent(mp);
        
        // get other relevant data for this integration point
        vec3d gcnt[3];
        vec3d ghr, ghs, ghrr, ghrs, ghss;
        vec3d tn, trn, tsn;
        double H, Hr, Hs, L, Lr, Ls;
        ContraBaseVectors(el, n, gcnt);
        MidShellSurface(el, n, L, Lr, Ls, H, Hr, Hs, ghr, ghs, ghrr, ghrs, ghss, tn, trn, tsn);
        double eta = el.gt(n);
        double gmag = (ghr ^ ghs).norm();
        mat3ds Imttg = (mat3dd(1) - dyad(tn))/gmag;
        mat3ds Imttgr = (dyads(tn, trn) + Imttg*(tn*((ghrr ^ ghs) + (ghr ^ ghrs))))/(-gmag);
        mat3ds Imttgs = (dyads(tn, tsn) + Imttg*(tn*((ghrs ^ ghs) + (ghr ^ ghss))))/(-gmag);
        mat3d Aghr, Aghs, Aghrr, Aghrs, Aghss;
        Aghr.skew(ghr); Aghs.skew(ghs);
        Aghrr.skew(ghrr); Aghrs.skew(ghrs); Aghss.skew(ghss);
        
        M = el.H(n); Mr = el.Hr(n); Ms = el.Hs(n);
        Mrr = el.Hrr(n); Mrs = el.Hrs(n); Mss = el.Hss(n);
        P = el.P(n); Pr = el.Pr(n); Ps = el.Ps(n);
        Prr = el.Prr(n); Prs = el.Prs(n); Pss = el.Pss(n);
        Q = el.Q(n); Qr = el.Qr(n); Qs = el.Qs(n);
        Qrr = el.Qrr(n); Qrs = el.Qrs(n); Qss = el.Qss(n);
        
        mat3ds Khr[4], Khs[4];
        mat3ds Khrr[4], Khrs[4], Khss[4];
        mat3d Qhr[4], Qhs[4];
        mat3d Qhrr[4], Qhrs[4], Qhss[4];
        mat3d Kr[4], Ks[4], Kt[4];
        mat3d dQr[4], dQs[4], dQt[4];
        
        mat3dd I(1);
        
        for (i=0; i<neln; ++i)
        {
            mat3ds pr, ps, prr, prs, pss;
            pr.zero(); ps.zero(); prr.zero(); prs.zero(); pss.zero();
            for (j=0; j<neln; ++j)
            {
                mat3ds Imtt = I - dyad(t[j]);
                pr += Imtt*(Pr[j]*fr[j][i] + Qr[j]*fs[j][i]);
                ps += Imtt*(Ps[j]*fr[j][i] + Qs[j]*fs[j][i]);
                prr += Imtt*(Prr[j]*fr[j][i] + Qrr[j]*fs[j][i]);
                prs += Imtt*(Prs[j]*fr[j][i] + Qrs[j]*fs[j][i]);
                pss += Imtt*(Pss[j]*fr[j][i] + Qss[j]*fs[j][i]);
            }
            
            Khr[i] = mat3dd(Mr[i]) + pr;
            Khs[i] = mat3dd(Ms[i]) + ps;
            Khrr[i] = mat3dd(Mrr[i]) + prr;
            Khrs[i] = mat3dd(Mrs[i]) + prs;
            Khss[i] = mat3dd(Mss[i]) + pss;
            Qhr[i] = Gr[i]*Pr[i] + Gs[i]*Qr[i];
            Qhs[i] = Gr[i]*Ps[i] + Gs[i]*Qs[i];
            Qhrr[i] = Gr[i]*Prr[i] + Gs[i]*Qrr[i];
            Qhrs[i] = Gr[i]*Prs[i] + Gs[i]*Qrs[i];
            Qhss[i] = Gr[i]*Pss[i] + Gs[i]*Qss[i];
            
            mat3d AK = Imttg*(Aghr*Khs[i] - Aghs*Khr[i]);
            Kr[i] = Khr[i] + AK*((L*Hr + Lr*H)*eta/2)
            + Imttgr*(Aghr*Khs[i] - Aghs*Khr[i])*(L*H*eta/2)
            + Imttg*(Aghrr*Khs[i] - Aghrs*Khr[i] + Aghr*Khrs[i] - Aghs*Khrr[i])*(L*H*eta/2);
            Ks[i] = Khs[i] + AK*((L*Hs + Ls*H)*eta/2)
            + Imttgs*(Aghr*Khs[i] - Aghs*Khr[i])*(L*H*eta/2)
            + Imttg*(Aghrs*Khs[i] - Aghss*Khr[i] + Aghr*Khss[i] - Aghs*Khrs[i])*(L*H*eta/2);
            Kt[i] = AK*(L*H/2);
            
            mat3d AQ = (tn & t[i])*M[i] + Imttg*(Aghr*Qhs[i] - Aghs*Qhr[i])*L;
            dQr[i] = Qhr[i] + AQ*(Hr*eta/2)
            + (((tn*Mr[i] + trn*M[i]) & t[i]) + (Imttg*Lr + Imttgr*L)*(Aghr*Qhs[i] - Aghs*Qhr[i]))*(H*eta/2)
            + Imttg*(Aghrr*Qhs[i] - Aghrs*Qhr[i] + Aghr*Qhrs[i] - Aghs*Qhrr[i])*(L*H*eta/2);
            dQs[i] = Qhs[i] + AQ*(Hs*eta/2)
            + (((tn*Ms[i] + tsn*M[i]) & t[i]) + (Imttg*Ls + Imttgs*L)*(Aghr*Qhs[i] - Aghs*Qhr[i]))*(H*eta/2)
            + Imttg*(Aghrs*Qhs[i] - Aghss*Qhr[i] + Aghr*Qhss[i] - Aghs*Qhrs[i])*(L*H*eta/2);
            dQt[i] = AQ*(H/2);
        }
        
        // ------------ constitutive component --------------
        
        // setup the material point
        // NOTE: deformation gradient has already been calculated in stress routine
        
        // we only calculate the upper triangular part
        // since ke is symmetric. The other part is
        // determined below using this symmetry.
        mat3d Kuu, Kud, Kdu, Kdd;
        for (i=0, i6=0; i<neln; ++i, i6 += 6)
        {
            for (j=0, j6 = 0; j<neln; ++j, j6 += 6)
            {
                Kuu
                = Kr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[0])*Kr[j]
                + Kr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[1])*Ks[j]
                + Kr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[2])*Kt[j]
                + Ks[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[0])*Kr[j]
                + Ks[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[1])*Ks[j]
                + Ks[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[2])*Kt[j]
                + Kt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[0])*Kr[j]
                + Kt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[1])*Ks[j]
                + Kt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[2])*Kt[j];
                
                Kud
                = Kr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[0])*dQr[j]
                + Kr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[1])*dQs[j]
                + Kr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[2])*dQt[j]
                + Ks[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[0])*dQr[j]
                + Ks[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[1])*dQs[j]
                + Ks[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[2])*dQt[j]
                + Kt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[0])*dQr[j]
                + Kt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[1])*dQs[j]
                + Kt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[2])*dQt[j];
                
                Kdu
                = dQr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[0])*Kr[j]
                + dQr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[1])*Ks[j]
                + dQr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[2])*Kt[j]
                + dQs[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[0])*Kr[j]
                + dQs[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[1])*Ks[j]
                + dQs[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[2])*Kt[j]
                + dQt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[0])*Kr[j]
                + dQt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[1])*Ks[j]
                + dQt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[2])*Kt[j];
                
                Kdd
                = dQr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[0])*dQr[j]
                + dQr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[1])*dQs[j]
                + dQr[i].transpose()*vdotTdotv(gcnt[0], C, gcnt[2])*dQt[j]
                + dQs[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[0])*dQr[j]
                + dQs[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[1])*dQs[j]
                + dQs[i].transpose()*vdotTdotv(gcnt[1], C, gcnt[2])*dQt[j]
                + dQt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[0])*dQr[j]
                + dQt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[1])*dQs[j]
                + dQt[i].transpose()*vdotTdotv(gcnt[2], C, gcnt[2])*dQt[j];
                
                Kuu *= detJ; Kud *= detJ; Kdu *= detJ; Kdd *= detJ;
                
                ke[i6  ][j6  ] += Kuu(0,0); ke[i6  ][j6+1] += Kuu(0,1); ke[i6  ][j6+2] += Kuu(0,2);
                ke[i6+1][j6  ] += Kuu(1,0); ke[i6+1][j6+1] += Kuu(1,1); ke[i6+1][j6+2] += Kuu(1,2);
                ke[i6+2][j6  ] += Kuu(2,0); ke[i6+2][j6+1] += Kuu(2,1); ke[i6+2][j6+2] += Kuu(2,2);
                
                ke[i6  ][j6+3] += Kud(0,0); ke[i6  ][j6+4] += Kud(0,1); ke[i6  ][j6+5] += Kud(0,2);
                ke[i6+1][j6+3] += Kud(1,0); ke[i6+1][j6+4] += Kud(1,1); ke[i6+1][j6+5] += Kud(1,2);
                ke[i6+2][j6+3] += Kud(2,0); ke[i6+2][j6+4] += Kud(2,1); ke[i6+2][j6+5] += Kud(2,2);
                
                ke[i6+3][j6  ] += Kdu(0,0); ke[i6+3][j6+1] += Kdu(0,1); ke[i6+3][j6+2] += Kdu(0,2);
                ke[i6+4][j6  ] += Kdu(1,0); ke[i6+4][j6+1] += Kdu(1,1); ke[i6+4][j6+2] += Kdu(1,2);
                ke[i6+5][j6  ] += Kdu(2,0); ke[i6+5][j6+1] += Kdu(2,1); ke[i6+5][j6+2] += Kdu(2,2);
                
                ke[i6+3][j6+3] += Kdd(0,0); ke[i6+3][j6+4] += Kdd(0,1); ke[i6+3][j6+5] += Kdd(0,2);
                ke[i6+4][j6+3] += Kdd(1,0); ke[i6+4][j6+4] += Kdd(1,1); ke[i6+4][j6+5] += Kdd(1,2);
                ke[i6+5][j6+3] += Kdd(2,0); ke[i6+5][j6+4] += Kdd(2,1); ke[i6+5][j6+5] += Kdd(2,2);
            }
        }
        
        // ------------ initial stress component --------------
        
        for (i=0, i6=0; i<neln; ++i, i6 += 6)
        {
            for (j=0, j6 = 0; j<neln; ++j, j6 += 6)
            {
                Kuu
                = Kr[i].transpose()*Kr[j]*(gcnt[0]*(T*gcnt[0]))
                + Kr[i].transpose()*Ks[j]*(gcnt[0]*(T*gcnt[1]))
                + Kr[i].transpose()*Kt[j]*(gcnt[0]*(T*gcnt[2]))
                + Ks[i].transpose()*Kr[j]*(gcnt[1]*(T*gcnt[0]))
                + Ks[i].transpose()*Ks[j]*(gcnt[1]*(T*gcnt[1]))
                + Ks[i].transpose()*Kt[j]*(gcnt[1]*(T*gcnt[2]))
                + Kt[i].transpose()*Kr[j]*(gcnt[2]*(T*gcnt[0]))
                + Kt[i].transpose()*Ks[j]*(gcnt[2]*(T*gcnt[1]))
                + Kt[i].transpose()*Kt[j]*(gcnt[2]*(T*gcnt[2]));
                
                Kud
                = Kr[i].transpose()*dQr[j]*(gcnt[0]*(T*gcnt[0]))
                + Kr[i].transpose()*dQs[j]*(gcnt[0]*(T*gcnt[1]))
                + Kr[i].transpose()*dQt[j]*(gcnt[0]*(T*gcnt[2]))
                + Ks[i].transpose()*dQr[j]*(gcnt[1]*(T*gcnt[0]))
                + Ks[i].transpose()*dQs[j]*(gcnt[1]*(T*gcnt[1]))
                + Ks[i].transpose()*dQt[j]*(gcnt[1]*(T*gcnt[2]))
                + Kt[i].transpose()*dQr[j]*(gcnt[2]*(T*gcnt[0]))
                + Kt[i].transpose()*dQs[j]*(gcnt[2]*(T*gcnt[1]))
                + Kt[i].transpose()*dQt[j]*(gcnt[2]*(T*gcnt[2]));
                
                Kdu
                = dQr[i].transpose()*Kr[j]*(gcnt[0]*(T*gcnt[0]))
                + dQr[i].transpose()*Ks[j]*(gcnt[0]*(T*gcnt[1]))
                + dQr[i].transpose()*Kt[j]*(gcnt[0]*(T*gcnt[2]))
                + dQs[i].transpose()*Kr[j]*(gcnt[1]*(T*gcnt[0]))
                + dQs[i].transpose()*Ks[j]*(gcnt[1]*(T*gcnt[1]))
                + dQs[i].transpose()*Kt[j]*(gcnt[1]*(T*gcnt[2]))
                + dQt[i].transpose()*Kr[j]*(gcnt[2]*(T*gcnt[0]))
                + dQt[i].transpose()*Ks[j]*(gcnt[2]*(T*gcnt[1]))
                + dQt[i].transpose()*Kt[j]*(gcnt[2]*(T*gcnt[2]));
                
                Kdd
                = dQr[i].transpose()*dQr[j]*(gcnt[0]*(T*gcnt[0]))
                + dQr[i].transpose()*dQs[j]*(gcnt[0]*(T*gcnt[1]))
                + dQr[i].transpose()*dQt[j]*(gcnt[0]*(T*gcnt[2]))
                + dQs[i].transpose()*dQr[j]*(gcnt[1]*(T*gcnt[0]))
                + dQs[i].transpose()*dQs[j]*(gcnt[1]*(T*gcnt[1]))
                + dQs[i].transpose()*dQt[j]*(gcnt[1]*(T*gcnt[2]))
                + dQt[i].transpose()*dQr[j]*(gcnt[2]*(T*gcnt[0]))
                + dQt[i].transpose()*dQs[j]*(gcnt[2]*(T*gcnt[1]))
                + dQt[i].transpose()*dQt[j]*(gcnt[2]*(T*gcnt[2]));
                
                Kuu *= detJ; Kud *= detJ; Kdu *= detJ; Kdd *= detJ;
                
                ke[i6  ][j6  ] += Kuu(0,0); ke[i6  ][j6+1] += Kuu(0,1); ke[i6  ][j6+2] += Kuu(0,2);
                ke[i6+1][j6  ] += Kuu(1,0); ke[i6+1][j6+1] += Kuu(1,1); ke[i6+1][j6+2] += Kuu(1,2);
                ke[i6+2][j6  ] += Kuu(2,0); ke[i6+2][j6+1] += Kuu(2,1); ke[i6+2][j6+2] += Kuu(2,2);
                
                ke[i6  ][j6+3] += Kud(0,0); ke[i6  ][j6+4] += Kud(0,1); ke[i6  ][j6+5] += Kud(0,2);
                ke[i6+1][j6+3] += Kud(1,0); ke[i6+1][j6+4] += Kud(1,1); ke[i6+1][j6+5] += Kud(1,2);
                ke[i6+2][j6+3] += Kud(2,0); ke[i6+2][j6+4] += Kud(2,1); ke[i6+2][j6+5] += Kud(2,2);
                
                ke[i6+3][j6  ] += Kdu(0,0); ke[i6+3][j6+1] += Kdu(0,1); ke[i6+3][j6+2] += Kdu(0,2);
                ke[i6+4][j6  ] += Kdu(1,0); ke[i6+4][j6+1] += Kdu(1,1); ke[i6+4][j6+2] += Kdu(1,2);
                ke[i6+5][j6  ] += Kdu(2,0); ke[i6+5][j6+1] += Kdu(2,1); ke[i6+5][j6+2] += Kdu(2,2);
                
                ke[i6+3][j6+3] += Kdd(0,0); ke[i6+3][j6+4] += Kdd(0,1); ke[i6+3][j6+5] += Kdd(0,2);
                ke[i6+4][j6+3] += Kdd(1,0); ke[i6+4][j6+4] += Kdd(1,1); ke[i6+4][j6+5] += Kdd(1,2);
                ke[i6+5][j6+3] += Kdd(2,0); ke[i6+5][j6+4] += Kdd(2,1); ke[i6+5][j6+5] += Kdd(2,2);
            }
        }
        
    } // end loop over gauss-points
}


//-----------------------------------------------------------------------------
//! Calculates body forces for shells

void FEElasticFergusonShellDomain::ElementBodyForce(FEModel& fem, FEFergusonShellElement& el, vector<double>& fe)
{
    int NF = fem.BodyLoads();
    for (int nf = 0; nf < NF; ++nf)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(nf));
        if (pbf)
        {
            double dens = m_pMat->Density();
            
            // calculate the average thickness
            double* h0 = &el.m_h0[0], gt, za;
            
            // integration weights
            double* gw = el.GaussWeights();
            double *Hn, detJ;
            
            // loop over integration points
            int nint = el.GaussPoints();
            int neln = el.Nodes();
            
            // nodal coordinates
            vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
            for (int i=0; i<neln; ++i)
            {
                r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
                rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
            }
            
            for (int n=0; n<nint; ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
                pt.m_r0 = el.Evaluate(r0, n);
                pt.m_rt = el.Evaluate(rt, n);
                
                detJ = detJ0(el, n)*gw[n];
                Hn  = el.H(n);
                gt = el.gt(n);
                
                // get the force
                vec3d f = pbf->force(mp);
                
                for (int i=0; i<neln; ++i)
                {
                    za = 0.5*gt*h0[i];
                    
                    fe[6*i  ] -= Hn[i]*f.x*dens*detJ;
                    fe[6*i+1] -= Hn[i]*f.y*dens*detJ;
                    fe[6*i+2] -= Hn[i]*f.z*dens*detJ;
                    
                    fe[6*i+3] -= za*Hn[i]*dens*f.x*detJ;
                    fe[6*i+4] -= za*Hn[i]*dens*f.y*detJ;
                    fe[6*i+5] -= za*Hn[i]*dens*f.z*detJ;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEElasticFergusonShellDomain::Update()
{
    FEMesh& mesh = *GetMesh();
    vec3d r0[9], rt[9];
    
    int n;
    for (int i=0; i<(int) m_Elem.size(); ++i)
    {
        // get the solid element
        FEFergusonShellElement& el = m_Elem[i];
        
        // get the number of integration points
        int nint = el.GaussPoints();
        
        // number of nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int j=0; j<neln; ++j)
        {
            r0[j] = mesh.Node(el.m_node[j]).m_r0;
            rt[j] = mesh.Node(el.m_node[j]).m_rt;
        }
        
        // get the integration weights
        double* gw = el.GaussWeights();
        
        // loop over the integration points and calculate
        // the stress at the integration point
        for (n=0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *(el.GetMaterialPoint(n));
            FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
            
            // material point coordinates
            // TODO: I'm not entirly happy with this solution
            //		 since the material point coordinates are used by most materials.
            pt.m_r0 = el.Evaluate(r0, n);
            pt.m_rt = el.Evaluate(rt, n);
            
            // get the deformation gradient and determinant
            pt.m_J = defgrad(el, pt.m_F, n);
            
            // calculate the stress at this material point
            pt.m_s = m_pMat->Stress(mp);
        }
    }
}


//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure
//! Note that for the shell elements the lm order is different compared
//! to the solid element ordering. This is because for shell elements the
//! nodes have six degrees of freedom each, where for solids they only
//! have 3 dofs.
void FEElasticFergusonShellDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
    int N = el.Nodes();
    lm.resize(N*9);
    for (int i=0; i<N; ++i)
    {
        FENode& node = m_pMesh->Node(el.m_node[i]);
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[6*i  ] = id[m_dofX];
        lm[6*i+1] = id[m_dofY];
        lm[6*i+2] = id[m_dofZ];
        
        // next the rotational dofs
        lm[6*i+3] = id[m_dofU];
        lm[6*i+4] = id[m_dofV];
        lm[6*i+5] = id[m_dofW];
        
        // rigid rotational dofs
        lm[6*N + 3*i  ] = id[m_dofRU];
        lm[6*N + 3*i+1] = id[m_dofRV];
        lm[6*N + 3*i+2] = id[m_dofRW];
    }
}
