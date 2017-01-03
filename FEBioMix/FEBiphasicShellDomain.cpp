//
//  FEBiphasicShellDomain.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 12/1/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FEBiphasicShellDomain.h"
#include "FECore/FEMesh.h"
#include "FECore/log.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
FEBiphasicShellDomain::FEBiphasicShellDomain(FEModel* pfem) : FESSIShellDomain(&pfem->GetMesh()), FEBiphasicDomain(pfem)
{
    m_dofU = pfem->GetDOFIndex("u");
    m_dofV = pfem->GetDOFIndex("v");
    m_dofW = pfem->GetDOFIndex("w");
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::SetMaterial(FEMaterial* pmat)
{
    m_pMat = dynamic_cast<FEBiphasic*>(pmat);
    assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEBiphasicShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    // initialize base class
    FESSIShellDomain::PreSolveUpdate(timeInfo);

    const int NE = FEElement::MAX_NODES;
    vec3d x0[NE], xt[NE], r0, rt;
    double pn[NE], qn[NE], p, q;
    FEMesh& m = *GetMesh();
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FEShellElement& el = m_Elem[i];
        int neln = el.Nodes();
        for (int i=0; i<neln; ++i)
        {
            x0[i] = m.Node(el.m_node[i]).m_r0;
            xt[i] = m.Node(el.m_node[i]).m_rt;
            pn[i] = m.Node(el.m_node[i]).get(m_dofP);
            qn[i] = m.Node(el.m_node[i]).get(m_dofQ);
        }
        
        int n = el.GaussPoints();
        for (int j=0; j<n; ++j)
        {
            r0 = el.Evaluate(x0, j);
            rt = el.Evaluate(xt, j);
            p = el.Evaluate(pn, j);
            q = el.Evaluate(qn, j);
            
            FEMaterialPoint& mp = *el.GetMaterialPoint(j);
            FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
            FEBiphasicMaterialPoint& pb = *mp.ExtractData<FEBiphasicMaterialPoint>();
            pt.m_r0 = r0;
            pt.m_rt = rt;
            
            pt.m_J = defgrad(el, pt.m_F, j);
            
            pb.m_Jp = pt.m_J;
            
            pb.m_p = p;
            pb.m_gradp = gradient(el, pn, qn, j);
            pb.m_phi0p = pb.m_phi0;
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
bool FEBiphasicShellDomain::Initialize()
{
    // initialize base class
    FESSIShellDomain::Initialize();
    
    // error flag (set true on error)
    bool bmerr = false;
    
    // initialize local coordinate systems (can I do this elsewhere?)
    FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FEShellElement& el = m_Elem[i];
        for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
    }
    
    // check for initially inverted shells
    for (int i=0; i<Elements(); ++i)
    {
        FEShellElement& el = Element(i);
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
void FEBiphasicShellDomain::Activate()
{
    for (int i=0; i<Nodes(); ++i)
    {
        FENode& node = Node(i);
        if (node.HasFlags(FENode::EXCLUDE) == false)
        {
            if (node.m_rid < 0)
            {
                node.m_ID[m_dofX] = DOF_ACTIVE;
                node.m_ID[m_dofY] = DOF_ACTIVE;
                node.m_ID[m_dofZ] = DOF_ACTIVE;
                
                if (node.HasFlags(FENode::SHELL))
                {
                    node.m_ID[m_dofU] = DOF_ACTIVE;
                    node.m_ID[m_dofV] = DOF_ACTIVE;
                    node.m_ID[m_dofW] = DOF_ACTIVE;
                }
            }
            
            node.m_ID[m_dofP] = DOF_ACTIVE;
            
            if (node.HasFlags(FENode::SHELL))
                node.m_ID[m_dofQ] = DOF_ACTIVE;
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data.
void FEBiphasicShellDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
    int N = el.Nodes();
    lm.resize(N*11);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = m_pMesh->Node(n);
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[8*i  ] = id[m_dofX];
        lm[8*i+1] = id[m_dofY];
        lm[8*i+2] = id[m_dofZ];
        
        // next the rotational dofs
        lm[8*i+3] = id[m_dofU];
        lm[8*i+4] = id[m_dofV];
        lm[8*i+5] = id[m_dofW];
        
        // now the pressure dofs
        lm[8*i+6] = id[m_dofP];
        lm[8*i+7] = id[m_dofQ];
        
        // rigid rotational dofs
        // TODO: Do I really need this?
        lm[8*N + 3*i  ] = id[m_dofRU];
        lm[8*N + 3*i+1] = id[m_dofRV];
        lm[8*N + 3*i+2] = id[m_dofRW];
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::Reset()
{
    // reset base class data
    FESSIShellDomain::Reset();
    
    // get the biphasic material
    FEBiphasic* pmb = m_pMat;
    
    // initialize all element data
    for (int i=0; i<(int) m_Elem.size(); ++i)
    {
        // get the solid element
        FEShellElement& el = m_Elem[i];
        
        // get the number of integration points
        int nint = el.GaussPoints();
        
        // loop over the integration points
        for (int n=0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            
            // initialize referential solid volume fraction
            pt.m_phi0 = pmb->m_phi0;
        }
    }
}


//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::InternalForces(FEGlobalVector& R)
{
    int NE = m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 8*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        //#pragma omp critical
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for shell elements

void FEBiphasicShellDomain::ElementInternalForce(FEShellElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix determinant
    double detJt;
    
    vec3d gradM, gradMu, gradMd;
    double Mu, Md;
    mat3ds s;
    
    const double* Mr, *Ms, *M;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    double*	gw = el.GaussWeights();
    double eta;
    
    vec3d gcnt[3];
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate the jacobian
        detJt = detJ(el, n);
        
        detJt *= gw[n];
        
        // get the stress vector for this integration point
        s = pt.m_s;
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        // next we get the determinant
        double Jp = bpt.m_Jp;
        double J = pt.m_J;
        
        // and then finally
        double divv = ((J-Jp)/dt)/J;
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        // get the solvent supply
        double phiwhat = m_pMat->SolventSupply(mp);
        
        for (i=0; i<neln; ++i)
        {
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMd = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            vec3d fu = s*gradMu;
            vec3d fd = s*gradMd;
            Mu = (1+eta)/2*M[i];
            Md = (1-eta)/2*M[i];
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[8*i  ] -= fu.x*detJt;
            fe[8*i+1] -= fu.y*detJt;
            fe[8*i+2] -= fu.z*detJt;
            
            fe[8*i+3] -= fd.x*detJt;
            fe[8*i+4] -= fd.y*detJt;
            fe[8*i+5] -= fd.z*detJt;
            
            fe[8*i+6] -= dt*(w*gradMu + (phiwhat - divv)*Mu)*detJt;
            fe[8*i+7] -= dt*(w*gradMd + (phiwhat - divv)*Md)*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::InternalForcesSS(FEGlobalVector& R)
{
    int NE = m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 8*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForceSS(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        //#pragma omp critical
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for shell elements
//! stead-state analysis

void FEBiphasicShellDomain::ElementInternalForceSS(FEShellElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix determinant
    double detJt;
    
    vec3d gradM, gradMu, gradMd;
    double Mu, Md;
    mat3ds s;
    
    const double* Mr, *Ms, *M;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    double*	gw = el.GaussWeights();
    double eta;
    
    vec3d gcnt[3];
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate the jacobian
        detJt = detJ(el, n);
        
        detJt *= gw[n];
        
        // get the stress vector for this integration point
        s = pt.m_s;
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        // get the solvent supply
        double phiwhat = m_pMat->SolventSupply(mp);
        
        for (i=0; i<neln; ++i)
        {
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMd = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            vec3d fu = s*gradMu;
            vec3d fd = s*gradMd;
            Mu = (1+eta)/2*M[i];
            Md = (1-eta)/2*M[i];
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[8*i  ] -= fu.x*detJt;
            fe[8*i+1] -= fu.y*detJt;
            fe[8*i+2] -= fu.z*detJt;
            
            fe[8*i+3] -= fd.x*detJt;
            fe[8*i+4] -= fd.y*detJt;
            fe[8*i+5] -= fd.z*detJt;
            
            fe[8*i+6] -= dt*(w*gradMu + phiwhat*Mu)*detJt;
            fe[8*i+7] -= dt*(w*gradMd + phiwhat*Md)*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::StiffnessMatrix(FESolver* psolver, bool bsymm)
{
    // repeat over all solid elements
    int NE = m_Elem.size();
    
#pragma omp parallel for shared(NE)
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FEShellElement& el = m_Elem[iel];
        
        // allocate stiffness matrix
        int neln = el.Nodes();
        int ndof = neln*8;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementBiphasicStiffness(el, ke, bsymm);
        
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::StiffnessMatrixSS(FESolver* psolver, bool bsymm)
{
    // repeat over all solid elements
    int NE = m_Elem.size();
    
#pragma omp parallel for shared(NE)
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FEShellElement& el = m_Elem[iel];
        
        // allocate stiffness matrix
        int neln = el.Nodes();
        int ndof = neln*8;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementBiphasicStiffnessSS(el, ke, bsymm);
        
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FEBiphasicShellDomain::ElementBiphasicStiffness(FEShellElement& el, matrix& ke, bool bsymm)
{
    int i, j, n;
    
    // jacobian matrix determinant
    double detJt;
    
    mat3ds s;
    
    const double* Mr, *Ms, *M;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    vector<vec3d> gradMu(neln), gradMd(neln);
    vector<double> Mu(neln), Md(neln);
    vec3d gradM;
    
    double*	gw = el.GaussWeights();
    double eta;
    
    vec3d gcnt[3];
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate the jacobian
        detJt = detJ(el, n);
        
        detJt *= gw[n];
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        for (i=0; i<neln; ++i)
        {
            // calculate global gradient of shape functions
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu[i] = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMd[i] = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            Mu[i] = (1+eta)/2*M[i];
            Md[i] = (1-eta)/2*M[i];
        }
        
        // get elasticity tensor
        tens4ds c = m_pMat->Tangent(pt);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = pt.m_gradp;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        // evaluate the solvent supply and its derivatives
        double phiwhat = 0;
        mat3ds Phie; Phie.zero();
        double Phip = 0;
        if (m_pMat->GetSolventSupply()) {
            phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
            Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
            Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
        }
        
        // Miscellaneous constants
        mat3dd I(1);
        
        // Kuu matrix
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                mat3d Kuu = (mat3dd(gradMu[i]*(s*gradMu[j])) + vdotTdotv(gradMu[i], c, gradMu[j]))*detJt;
                mat3d Kud = (mat3dd(gradMu[i]*(s*gradMd[j])) + vdotTdotv(gradMu[i], c, gradMd[j]))*detJt;
                mat3d Kdu = (mat3dd(gradMd[i]*(s*gradMu[j])) + vdotTdotv(gradMd[i], c, gradMu[j]))*detJt;
                mat3d Kdd = (mat3dd(gradMd[i]*(s*gradMd[j])) + vdotTdotv(gradMd[i], c, gradMd[j]))*detJt;
                
                ke[8*i  ][8*j  ] += Kuu[0][0]; ke[8*i  ][8*j+1] += Kuu[0][1]; ke[8*i  ][8*j+2] += Kuu[0][2];
                ke[8*i+1][8*j  ] += Kuu[1][0]; ke[8*i+1][8*j+1] += Kuu[1][1]; ke[8*i+1][8*j+2] += Kuu[1][2];
                ke[8*i+2][8*j  ] += Kuu[2][0]; ke[8*i+2][8*j+1] += Kuu[2][1]; ke[8*i+2][8*j+2] += Kuu[2][2];
                
                ke[8*i  ][8*j+3] += Kud[0][0]; ke[8*i  ][8*j+4] += Kud[0][1]; ke[8*i  ][8*j+5] += Kud[0][2];
                ke[8*i+1][8*j+3] += Kud[1][0]; ke[8*i+1][8*j+4] += Kud[1][1]; ke[8*i+1][8*j+5] += Kud[1][2];
                ke[8*i+2][8*j+3] += Kud[2][0]; ke[8*i+2][8*j+4] += Kud[2][1]; ke[8*i+2][8*j+5] += Kud[2][2];
                
                ke[8*i+3][8*j  ] += Kdu[0][0]; ke[8*i+3][8*j+1] += Kdu[0][1]; ke[8*i+3][8*j+2] += Kdu[0][2];
                ke[8*i+4][8*j  ] += Kdu[1][0]; ke[8*i+4][8*j+1] += Kdu[1][1]; ke[8*i+4][8*j+2] += Kdu[1][2];
                ke[8*i+5][8*j  ] += Kdu[2][0]; ke[8*i+5][8*j+1] += Kdu[2][1]; ke[8*i+5][8*j+2] += Kdu[2][2];
                
                ke[8*i+3][8*j+3] += Kdd[0][0]; ke[8*i+3][8*j+4] += Kdd[0][1]; ke[8*i+3][8*j+5] += Kdd[0][2];
                ke[8*i+4][8*j+3] += Kdd[1][0]; ke[8*i+4][8*j+4] += Kdd[1][1]; ke[8*i+4][8*j+5] += Kdd[1][2];
                ke[8*i+5][8*j+3] += Kdd[2][0]; ke[8*i+5][8*j+4] += Kdd[2][1]; ke[8*i+5][8*j+5] += Kdd[2][2];
                
            }
        
        // kup matrix
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                vec3d kup = gradMu[i]*(Mu[j]*detJt);
                vec3d kuq = gradMu[i]*(Md[j]*detJt);
                vec3d kdp = gradMd[i]*(Mu[j]*detJt);
                vec3d kdq = gradMd[i]*(Md[j]*detJt);
                
                ke[8*i  ][8*j+6] -= kup.z;
                ke[8*i+1][8*j+6] -= kup.y;
                ke[8*i+2][8*j+6] -= kup.z;
                
                ke[8*i  ][8*j+7] -= kuq.x;
                ke[8*i+1][8*j+7] -= kuq.y;
                ke[8*i+2][8*j+7] -= kuq.z;
                
                ke[8*i+3][8*j+6] -= kdp.x;
                ke[8*i+4][8*j+6] -= kdp.y;
                ke[8*i+5][8*j+6] -= kdp.z;
                
                ke[8*i+3][8*j+7] -= kdq.x;
                ke[8*i+4][8*j+7] -= kdq.y;
                ke[8*i+5][8*j+7] -= kdq.z;
                
            }
        
        // kpu matrix
        mat3ds Q = mat3dd(1/dt-phiwhat) - Phie*ept.m_J;
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                vec3d kpu = (vdotTdotv(gradMu[i], dKdE, gradMu[j]).transpose()*gradp + Q*(gradMu[j]*Mu[i]))*(detJt*dt);
                vec3d kpd = (vdotTdotv(gradMu[i], dKdE, gradMd[j]).transpose()*gradp + Q*(gradMd[j]*Mu[i]))*(detJt*dt);
                vec3d kqu = (vdotTdotv(gradMd[i], dKdE, gradMu[j]).transpose()*gradp + Q*(gradMu[j]*Md[i]))*(detJt*dt);
                vec3d kqd = (vdotTdotv(gradMd[i], dKdE, gradMd[j]).transpose()*gradp + Q*(gradMd[j]*Md[i]))*(detJt*dt);
                
                ke[8*i+6][8*j  ] -= kpu.x; ke[8*i+6][8*j+1] -= kpu.y; ke[8*i+6][8*j+2] -= kpu.z;
                
                ke[8*i+6][8*j+3] -= kpd.x; ke[8*i+6][8*j+4] -= kpd.y; ke[8*i+6][8*j+5] -= kpd.z;
                
                ke[8*i+7][8*j  ] -= kqu.x; ke[8*i+7][8*j+1] -= kqu.y; ke[8*i+7][8*j+2] -= kqu.z;
                
                ke[8*i+7][8*j+3] -= kqd.x; ke[8*i+7][8*j+4] -= kqd.y; ke[8*i+7][8*j+5] -= kqd.z;
                
            }
        
        // kpp matrix
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                double kpp = (gradMu[i]*(K*gradMu[j]) - Phip*Mu[i]*Mu[j])*(detJt*dt);
                double kpq = (gradMu[i]*(K*gradMd[j]) - Phip*Mu[i]*Mu[j])*(detJt*dt);
                double kqp = (gradMd[i]*(K*gradMu[j]) - Phip*Md[i]*Mu[j])*(detJt*dt);
                double kqq = (gradMd[i]*(K*gradMd[j]) - Phip*Md[i]*Md[j])*(detJt*dt);
                
                ke[8*i+6][8*j+6] -= kpp;
                
                ke[8*i+6][8*j+7] -= kpq;
                
                ke[8*i+7][8*j+6] -= kqp;
                
                ke[8*i+7][8*j+7] -= kqq;
                
            }
        
    }
    return true;
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//! for the steady-state response (zero solid velocity)
//!
bool FEBiphasicShellDomain::ElementBiphasicStiffnessSS(FEShellElement& el, matrix& ke, bool bsymm)
{
    int i, j, n;
    
    // jacobian matrix determinant
    double detJt;
    
    mat3ds s;
    
    const double* Mr, *Ms, *M;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    vector<vec3d> gradMu(neln), gradMd(neln);
    vector<double> Mu(neln), Md(neln);
    vec3d gradM;
    
    double*	gw = el.GaussWeights();
    double eta;
    
    vec3d gcnt[3];
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate the jacobian
        detJt = detJ(el, n);
        
        detJt *= gw[n];
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        for (i=0; i<neln; ++i)
        {
            // calculate global gradient of shape functions
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu[i] = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMd[i] = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            Mu[i] = (1+eta)/2*M[i];
            Md[i] = (1-eta)/2*M[i];
        }
        
        // get elasticity tensor
        tens4ds c = m_pMat->Tangent(pt);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = pt.m_gradp;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        // evaluate the solvent supply and its derivatives
        double phiwhat = 0;
        mat3ds Phie; Phie.zero();
        double Phip = 0;
        if (m_pMat->GetSolventSupply()) {
            phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
            Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
            Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
        }
        
        // Miscellaneous constants
        mat3dd I(1);
        
        // Kuu matrix
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                mat3d Kuu = (mat3dd(gradMu[i]*(s*gradMu[j])) + vdotTdotv(gradMu[i], c, gradMu[j]))*detJt;
                mat3d Kud = (mat3dd(gradMu[i]*(s*gradMd[j])) + vdotTdotv(gradMu[i], c, gradMd[j]))*detJt;
                mat3d Kdu = (mat3dd(gradMd[i]*(s*gradMu[j])) + vdotTdotv(gradMd[i], c, gradMu[j]))*detJt;
                mat3d Kdd = (mat3dd(gradMd[i]*(s*gradMd[j])) + vdotTdotv(gradMd[i], c, gradMd[j]))*detJt;
                
                ke[8*i  ][8*j  ] += Kuu[0][0]; ke[8*i  ][8*j+1] += Kuu[0][1]; ke[8*i  ][8*j+2] += Kuu[0][2];
                ke[8*i+1][8*j  ] += Kuu[1][0]; ke[8*i+1][8*j+1] += Kuu[1][1]; ke[8*i+1][8*j+2] += Kuu[1][2];
                ke[8*i+2][8*j  ] += Kuu[2][0]; ke[8*i+2][8*j+1] += Kuu[2][1]; ke[8*i+2][8*j+2] += Kuu[2][2];
                
                ke[8*i  ][8*j+3] += Kud[0][0]; ke[8*i  ][8*j+4] += Kud[0][1]; ke[8*i  ][8*j+5] += Kud[0][2];
                ke[8*i+1][8*j+3] += Kud[1][0]; ke[8*i+1][8*j+4] += Kud[1][1]; ke[8*i+1][8*j+5] += Kud[1][2];
                ke[8*i+2][8*j+3] += Kud[2][0]; ke[8*i+2][8*j+4] += Kud[2][1]; ke[8*i+2][8*j+5] += Kud[2][2];
                
                ke[8*i+3][8*j  ] += Kdu[0][0]; ke[8*i+3][8*j+1] += Kdu[0][1]; ke[8*i+3][8*j+2] += Kdu[0][2];
                ke[8*i+4][8*j  ] += Kdu[1][0]; ke[8*i+4][8*j+1] += Kdu[1][1]; ke[8*i+4][8*j+2] += Kdu[1][2];
                ke[8*i+5][8*j  ] += Kdu[2][0]; ke[8*i+5][8*j+1] += Kdu[2][1]; ke[8*i+5][8*j+2] += Kdu[2][2];
                
                ke[8*i+3][8*j+3] += Kdd[0][0]; ke[8*i+3][8*j+4] += Kdd[0][1]; ke[8*i+3][8*j+5] += Kdd[0][2];
                ke[8*i+4][8*j+3] += Kdd[1][0]; ke[8*i+4][8*j+4] += Kdd[1][1]; ke[8*i+4][8*j+5] += Kdd[1][2];
                ke[8*i+5][8*j+3] += Kdd[2][0]; ke[8*i+5][8*j+4] += Kdd[2][1]; ke[8*i+5][8*j+5] += Kdd[2][2];
                
            }
        
        // kup matrix
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                vec3d kup = gradMu[i]*(Mu[j]*detJt);
                vec3d kuq = gradMu[i]*(Md[j]*detJt);
                vec3d kdp = gradMd[i]*(Mu[j]*detJt);
                vec3d kdq = gradMd[i]*(Md[j]*detJt);
                
                ke[8*i  ][8*j+6] -= kup.z;
                ke[8*i+1][8*j+6] -= kup.y;
                ke[8*i+2][8*j+6] -= kup.z;
                
                ke[8*i  ][8*j+7] -= kuq.x;
                ke[8*i+1][8*j+7] -= kuq.y;
                ke[8*i+2][8*j+7] -= kuq.z;
                
                ke[8*i+3][8*j+6] -= kdp.x;
                ke[8*i+4][8*j+6] -= kdp.y;
                ke[8*i+5][8*j+6] -= kdp.z;
                
                ke[8*i+3][8*j+7] -= kdq.x;
                ke[8*i+4][8*j+7] -= kdq.y;
                ke[8*i+5][8*j+7] -= kdq.z;
                
            }
        
        // kpu matrix
        mat3ds Q = mat3dd(-phiwhat) - Phie*ept.m_J;
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                vec3d kpu = (vdotTdotv(gradMu[i], dKdE, gradMu[j]).transpose()*gradp + Q*(gradMu[j]*Mu[i]))*(detJt*dt);
                vec3d kpd = (vdotTdotv(gradMu[i], dKdE, gradMd[j]).transpose()*gradp + Q*(gradMd[j]*Mu[i]))*(detJt*dt);
                vec3d kqu = (vdotTdotv(gradMd[i], dKdE, gradMu[j]).transpose()*gradp + Q*(gradMu[j]*Md[i]))*(detJt*dt);
                vec3d kqd = (vdotTdotv(gradMd[i], dKdE, gradMd[j]).transpose()*gradp + Q*(gradMd[j]*Md[i]))*(detJt*dt);
                
                ke[8*i+6][8*j  ] -= kpu.x; ke[8*i+6][8*j+1] -= kpu.y; ke[8*i+6][8*j+2] -= kpu.z;
                
                ke[8*i+6][8*j+3] -= kpd.x; ke[8*i+6][8*j+4] -= kpd.y; ke[8*i+6][8*j+5] -= kpd.z;
                
                ke[8*i+7][8*j  ] -= kqu.x; ke[8*i+7][8*j+1] -= kqu.y; ke[8*i+7][8*j+2] -= kqu.z;
                
                ke[8*i+7][8*j+3] -= kqd.x; ke[8*i+7][8*j+4] -= kqd.y; ke[8*i+7][8*j+5] -= kqd.z;
                
            }
        
        // kpp matrix
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                double kpp = (gradMu[i]*(K*gradMu[j]) - Phip*Mu[i]*Mu[j])*(detJt*dt);
                double kpq = (gradMu[i]*(K*gradMd[j]) - Phip*Mu[i]*Mu[j])*(detJt*dt);
                double kqp = (gradMd[i]*(K*gradMu[j]) - Phip*Md[i]*Mu[j])*(detJt*dt);
                double kqq = (gradMd[i]*(K*gradMd[j]) - Phip*Md[i]*Md[j])*(detJt*dt);
                
                ke[8*i+6][8*j+6] -= kpp;
                
                ke[8*i+6][8*j+7] -= kpq;
                
                ke[8*i+7][8*j+6] -= kqp;
                
                ke[8*i+7][8*j+7] -= kqq;
                
            }
        
    }
    return true;
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::Update(const FETimeInfo& tp)
{
    bool berr = false;
    int NE = (int) m_Elem.size();
#pragma omp parallel for shared(NE, berr)
    for (int i=0; i<NE; ++i)
    {
        try
        {
            UpdateElementStress(i);
        }
        catch (NegativeJacobian e)
        {
#pragma omp critical
            {
                berr = true;
                if (NegativeJacobian::m_boutput) e.print();
            }
        }
    }
    // if we encountered an error, we request a running restart
    if (berr)
    {
        if (NegativeJacobian::m_boutput == false) felog.printbox("ERROR", "Negative jacobian was detected.");
        throw DoRunningRestart();
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::UpdateElementStress(int iel)
{
    // get the solid element
    FEShellElement& el = m_Elem[iel];
    
    // get the number of integration points
    int nint = el.GaussPoints();
    
    // get the number of nodes
    int neln = el.Nodes();
    
    // get the nodal data
    FEMesh& mesh = *m_pMesh;
    vec3d r0[FEElement::MAX_NODES];
    vec3d rt[FEElement::MAX_NODES];
    double pn[FEElement::MAX_NODES];
    double qn[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j)
    {
        r0[j] = mesh.Node(el.m_node[j]).m_r0;
        rt[j] = mesh.Node(el.m_node[j]).m_rt;
        pn[j] = mesh.Node(el.m_node[j]).get(m_dofP);
        qn[j] = mesh.Node(el.m_node[j]).get(m_dofQ);
    }
    
    // loop over the integration points and calculate
    // the stress at the integration point
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        
        // material point coordinates
        // TODO: I'm not entirly happy with this solution
        //		 since the material point coordinates are used by most materials.
        pt.m_r0 = el.Evaluate(r0, n);
        pt.m_rt = el.Evaluate(rt, n);
        
        // get the deformation gradient and determinant
        pt.m_J = defgrad(el, pt.m_F, n);
        
        // poroelasticity data
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // evaluate fluid pressure at gauss-point
        ppt.m_p = evaluate(el, pn, qn, n);
        
        // calculate the gradient of p at gauss-point
        ppt.m_gradp = gradient(el, pn, qn, n);
        
        // for biphasic materials also update the fluid flux
        //		ppt.m_w = m_pMat->Flux(mp);
        ppt.m_w = FluidFlux(mp);
        ppt.m_pa = m_pMat->Pressure(mp);
        
        // calculate the stress at this material point
        pt.m_s = m_pMat->Stress(mp);
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 8*el.Nodes();
        fe.assign(ndof, 0);
        
        // apply body forces
        ElementBodyForce(BF, el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEBiphasicShellDomain::ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe)
{
    // get true solid and fluid densities
    double rhoTs = m_pMat->SolidDensity();
    double rhoTw = m_pMat->FluidDensity();
    
    // jacobian
    double detJt, eta;
    double *M;
    double* gw = el.GaussWeights();
    vec3d b, fu, fd;
    
    // number of nodes
    int neln = el.Nodes();
    
    // nodal coordinates
    vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
        rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
    }
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
        pt.m_r0 = el.Evaluate(r0, n);
        pt.m_rt = el.Evaluate(rt, n);
        
        detJt = detJ(el, n)*gw[n];
        
        // get the force
        b = BF.force(mp);
        
        eta = el.gt(n);
        
        // evaluate apparent solid and fluid densities and mixture density
        double phiw = m_pMat->Porosity(mp);
        double rhos = (1-phiw)*rhoTs;
        double rhow = phiw*rhoTw;
        double rho = rhos + rhow;
        
        M = el.H(n);
        
        for (int i=0; i<neln; ++i)
        {
            fu = b*(rho*(1+eta)/2*M[i]*detJt);
            fd = b*(rho*(1-eta)/2*M[i]*detJt);
            fe[8*i  ] -= fu.x;
            fe[8*i+1] -= fu.y;
            fe[8*i+2] -= fu.z;
            
            fe[8*i+3] -= fd.x;
            fe[8*i+4] -= fd.y;
            fe[8*i+5] -= fd.z;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::BodyForceStiffness(FESolver* psolver, FEBodyForce& bf)
{
    FEBiphasic* pmb = dynamic_cast<FEBiphasic*>(GetMaterial()); assert(pmb);
    
    // element stiffness matrix
    matrix ke;
    vector<int> lm;
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    for (int iel=0; iel<NE; ++iel)
    {
        FEShellElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int neln = el.Nodes();
        int ndof = 8*neln;
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementBodyForceStiffness(bf, el, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEBiphasicShellDomain::ElementBodyForceStiffness(FEBodyForce& BF, FEShellElement &el, matrix &ke)
{
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // get true solid and fluid densities
    double rhoTs = m_pMat->SolidDensity();
    double rhoTw = m_pMat->FluidDensity();
    
    // jacobian
    double detJt;
    const double* Mr, *Ms, *M;
    vector<vec3d> gradMu(neln), gradMd(neln);
    vector<double> Mu(neln), Md(neln);
    vec3d gradM;
    
    double* gw = el.GaussWeights();
    double eta;
    
    vec3d gcnt[3];
    
    vec3d b;
    mat3ds gradb;
    mat3d Kuu, Kud, Kdu, Kdd;
    vec3d kpu, kpd, kqu, kqd;
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        
        // get the body force
        b = BF.force(mp);
        
        // get the body force stiffness
        gradb = BF.stiffness(mp);
        
        // evaluate apparent solid and fluid densities and mixture density
        double phiw = m_pMat->Porosity(mp);
        double rhos = (1-phiw)*rhoTs;
        double rhow = phiw*rhoTw;
        double rho = rhos + rhow;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        // calculate the jacobian
        detJt = detJ(el, n);
        
        detJt *= gw[n];
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        for (int i=0; i<neln; ++i)
        {
            // calculate global gradient of shape functions
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu[i] = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMd[i] = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            Mu[i] = (1+eta)/2*M[i];
            Md[i] = (1-eta)/2*M[i];
        }
        
        for (int i=0; i<neln; ++i)
            for (int j=0; j<neln; ++j)
            {
                Kuu = (gradb*(Mu[j]*rho) + (b & gradMu[j])*rhoTw)*(Mu[i]*detJt);
                Kud = (gradb*(Md[j]*rho) + (b & gradMd[j])*rhoTw)*(Mu[i]*detJt);
                Kdu = (gradb*(Mu[j]*rho) + (b & gradMu[j])*rhoTw)*(Md[i]*detJt);
                Kdd = (gradb*(Md[j]*rho) + (b & gradMd[j])*rhoTw)*(Md[i]*detJt);
                
                ke[8*i  ][8*j  ] += Kuu(0,0); ke[8*i  ][8*j+1] += Kuu(0,1); ke[8*i  ][8*j+2] += Kuu(0,2);
                ke[8*i+1][8*j  ] += Kuu(1,0); ke[8*i+1][8*j+1] += Kuu(1,1); ke[8*i+1][8*j+2] += Kuu(1,2);
                ke[8*i+2][8*j  ] += Kuu(2,0); ke[8*i+2][8*j+1] += Kuu(2,1); ke[8*i+2][8*j+2] += Kuu(2,2);
                
                ke[8*i  ][8*j+3] += Kud(0,0); ke[8*i  ][8*j+4] += Kud(0,1); ke[8*i  ][8*j+5] += Kud(0,2);
                ke[8*i+1][8*j+3] += Kud(1,0); ke[8*i+1][8*j+4] += Kud(1,1); ke[8*i+1][8*j+5] += Kud(1,2);
                ke[8*i+2][8*j+3] += Kud(2,0); ke[8*i+2][8*j+4] += Kud(2,1); ke[8*i+2][8*j+5] += Kud(2,2);
                
                ke[8*i+3][8*j  ] += Kdu(0,0); ke[8*i+3][8*j+1] += Kdu(0,1); ke[8*i+3][8*j+2] += Kdu(0,2);
                ke[8*i+4][8*j  ] += Kdu(1,0); ke[8*i+4][8*j+1] += Kdu(1,1); ke[8*i+4][8*j+2] += Kdu(1,2);
                ke[8*i+5][8*j  ] += Kdu(2,0); ke[8*i+5][8*j+1] += Kdu(2,1); ke[8*i+5][8*j+2] += Kdu(2,2);
                
                ke[8*i+3][8*j+3] += Kdd(0,0); ke[8*i+3][8*j+4] += Kdd(0,1); ke[8*i+3][8*j+5] += Kdd(0,2);
                ke[8*i+4][8*j+3] += Kdd(1,0); ke[8*i+4][8*j+4] += Kdd(1,1); ke[8*i+4][8*j+5] += Kdd(1,2);
                ke[8*i+5][8*j+3] += Kdd(2,0); ke[8*i+5][8*j+4] += Kdd(2,1); ke[8*i+5][8*j+5] += Kdd(2,2);
                
                kpu = (vdotTdotv(gradMu[i], dKdE, gradMu[j]).transpose()*b
                       + ((b & gradMu[j]) + gradb*Mu[j])*K*gradMu[i])*(rhoTw*detJt*dt);
                kpd = (vdotTdotv(gradMu[i], dKdE, gradMd[j]).transpose()*b
                       + ((b & gradMd[j]) + gradb*Md[j])*K*gradMu[i])*(rhoTw*detJt*dt);
                kqu = (vdotTdotv(gradMd[i], dKdE, gradMu[j]).transpose()*b
                       + ((b & gradMu[j]) + gradb*Mu[j])*K*gradMd[i])*(rhoTw*detJt*dt);
                kqd = (vdotTdotv(gradMd[i], dKdE, gradMd[j]).transpose()*b
                       + ((b & gradMd[j]) + gradb*Md[j])*K*gradMd[i])*(rhoTw*detJt*dt);
                
                ke[8*i+6][8*j  ] -= kpu.x; ke[8*i+6][8*j+1] -= kpu.y; ke[8*i+6][8*j+2] -= kpu.z;
                ke[8*i+6][8*j+3] -= kpd.x; ke[8*i+6][8*j+4] -= kpd.y; ke[8*i+6][8*j+5] -= kpd.z;
                ke[8*i+7][8*j  ] -= kqu.x; ke[8*i+7][8*j+1] -= kqu.y; ke[8*i+7][8*j+2] -= kqu.z;
                ke[8*i+7][8*j+3] -= kqd.x; ke[8*i+7][8*j+4] -= kqd.y; ke[8*i+7][8*j+5] -= kqd.z;
            }
    }
    
}

//-----------------------------------------------------------------------------
vec3d FEBiphasicShellDomain::FluidFlux(FEMaterialPoint& mp)
{
    FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
    
    // pressure gradient
    vec3d gradp = ppt.m_gradp;
    
    // fluid flux w = -k*grad(p)
    mat3ds kt = m_pMat->Permeability(mp);
    
    vec3d w = -(kt*gradp);
    
    // get true fluid density
    double rhoTw = m_pMat->FluidDensity();
    
    // body force contribution
    FEModel& fem = *m_pMat->GetFEModel();
    int nbf = fem.BodyLoads();
    if (nbf) {
        vec3d b(0,0,0);
        for (int i=0; i<nbf; ++i)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(i));
            if (pbf->IsActive())
            {
                // negate b because body forces are defined with a negative sign in FEBio
                b -= pbf->force(mp);
            }
        }
        w += (kt*b)*(rhoTw);
    }
    
    // active momentum supply contribution
    FEActiveMomentumSupply* pAmom = m_pMat->GetActiveMomentumSupply();
    if (pAmom) {
        vec3d pw = pAmom->ActiveSupply(mp);
        w += kt*pw;
    }
    
    return w;
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FEBiphasicShellDomain::CoBaseVectors(FEShellElement& el, int n, vec3d g[3])
{
    int i;
    
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r[i] = ni.m_rt;
        D[i] = ni.m_d0 + ni.get_vec3d(m_dofX, m_dofY, m_dofZ) - ni.get_vec3d(m_dofU, m_dofV, m_dofW);
    }
    
    double eta = el.gt(n);
    
    double* Mr = el.Hr(n);
    double* Ms = el.Hs(n);
    double* M  = el.H(n);
    
    // initialize covariant basis vectors
    g[0] = g[1] = g[2] = vec3d(0,0,0);
    
    for (i=0; i<neln; ++i)
    {
        g[0] += (r[i] - D[i]*(1-eta)/2)*Mr[i];
        g[1] += (r[i] - D[i]*(1-eta)/2)*Ms[i];
        g[2] += D[i]*(M[i]/2);
    }
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FEBiphasicShellDomain::ContraBaseVectors(FEShellElement& el, int n, vec3d gcnt[3])
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
// jacobian with respect to current frame
double FEBiphasicShellDomain::detJ(FEShellElement& el, int n)
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    return J.det();
}

//-----------------------------------------------------------------------------
double FEBiphasicShellDomain::defgrad(FEShellElement& el, mat3d& F, int n)
{
    vec3d gcov[3], Gcnt[3];
    CoBaseVectors(el, n, gcov);
    ContraBaseVectors0(el, n, Gcnt);
    
    F = (gcov[0] & Gcnt[0]) + (gcov[1] & Gcnt[1]) + (gcov[2] & Gcnt[2]);
    double J = F.det();
    if (J <= 0) throw NegativeJacobian(el.GetID(), n, J, &el);
    
    return J;
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the current frame at
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
double FEBiphasicShellDomain::invjact(FEShellElement& el, double Ji[3][3], int n)
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
//! evaluate a scalar function over the shell
double FEBiphasicShellDomain::evaluate(FEShellElement& el, double* pn, double* dpn, int n)
{
    const double *M;
    double eta;
    vec3d gcnt[3];
    double Mu, Md;
    double p = 0;
    
    eta = el.gt(n);
    
    M  = el.H(n);
    
    ContraBaseVectors(el, n, gcnt);
    
    int neln = el.Nodes();
    
    for (int i=0; i<neln; ++i)
    {
        Mu = (1+eta)/2*M[i];
        Md = (1-eta)/2*M[i];
        p += Mu*pn[i] + Md*dpn[i];
    }
    
    return p;
}

//-----------------------------------------------------------------------------
//! calculate the gradient of a scalar function over the shell
vec3d FEBiphasicShellDomain::gradient(FEShellElement& el, double* pn, double* dpn, int n)
{
    const double* Mr, *Ms, *M;
    double eta;
    vec3d gcnt[3];
    vec3d gradM, gradMu, gradMd;
    vec3d gradp(0,0,0);
    
    eta = el.gt(n);
    
    Mr = el.Hr(n);
    Ms = el.Hs(n);
    M  = el.H(n);
    
    ContraBaseVectors(el, n, gcnt);
    
    int neln = el.Nodes();
    
    for (int i=0; i<neln; ++i)
    {
        gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
        gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
        gradMd = (gradM*(1-eta) - gcnt[2]*M[i])/2;
        gradp += gradMu*pn[i] + gradMd*dpn[i];
    }
    
    return gradp;
}
