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
#include "FEBiphasicShellDomain.h"
#include "FECore/FEMesh.h"
#include "FECore/log.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
FEBiphasicShellDomain::FEBiphasicShellDomain(FEModel* pfem) : FESSIShellDomain(pfem), FEBiphasicDomain(pfem), m_dof(pfem)
{
    m_dofSX = pfem->GetDOFIndex("sx");
    m_dofSY = pfem->GetDOFIndex("sy");
    m_dofSZ = pfem->GetDOFIndex("sz");
}

//-----------------------------------------------------------------------------
//! get the material (overridden from FEDomain)
FEMaterial* FEBiphasicShellDomain::GetMaterial()
{
	return m_pMat;
}

//-----------------------------------------------------------------------------
//! get the total dof
const FEDofList& FEBiphasicShellDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
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
    double pn[NE], qn[NE], p;
    FEMesh& m = *GetMesh();
    for (size_t iel=0; iel<m_Elem.size(); ++iel)
    {
        FEShellElement& el = m_Elem[iel];
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
            p = evaluate(el, pn, qn, j);
            
            FEMaterialPoint& mp = *el.GetMaterialPoint(j);
            FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
            FEBiphasicMaterialPoint& pb = *mp.ExtractData<FEBiphasicMaterialPoint>();
            mp.m_r0 = r0;
            mp.m_rt = rt;
            
            pt.m_J = defgrad(el, pt.m_F, j);
            
            pb.m_Jp = pt.m_J;
            
            pb.m_p = p;
            pb.m_gradp = gradient(el, pn, qn, j);
            pb.m_gradpp = pb.m_gradp;
            pb.m_phi0p = pb.m_phi0t;
            
            mp.Update(timeInfo);
        }
    }
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
                node.set_active(m_dofU[0]);
                node.set_active(m_dofU[1]);
                node.set_active(m_dofU[2]);
                
                if (node.HasFlags(FENode::SHELL))
                {
                    node.set_active(m_dofSU[0]);
                    node.set_active(m_dofSU[1]);
                    node.set_active(m_dofSU[2]);
                }
            }
            
            node.set_active(m_dofP);
            
            if (node.HasFlags(FENode::SHELL))
                node.set_active(m_dofQ);
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
        lm[8*i  ] = id[m_dofU[0]];
        lm[8*i+1] = id[m_dofU[1]];
        lm[8*i+2] = id[m_dofU[2]];
        
        // next the shell dofs
        lm[8*i+3] = id[m_dofSX];
        lm[8*i+4] = id[m_dofSY];
        lm[8*i+5] = id[m_dofSZ];
        
        // now the pressure dofs
        lm[8*i+6] = id[m_dofP];
        lm[8*i+7] = id[m_dofQ];
        
        // rigid rotational dofs
        // TODO: Do I really need this?
        lm[8*N + 3*i  ] = id[m_dofR[0]];
        lm[8*N + 3*i+1] = id[m_dofR[1]];
        lm[8*N + 3*i+2] = id[m_dofR[2]];
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
	ForEachMaterialPoint([=](FEMaterialPoint& mp) {
		FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());

		// initialize referential solid volume fraction
		pt.m_phi0 = pt.m_phi0t = pmb->m_phi0(mp);
	});
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::InternalForces(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
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
        R.Assemble(el.m_node, lm, fe, true);
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
    int NE = (int)m_Elem.size();
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
        R.Assemble(el.m_node, lm, fe, true);
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
void FEBiphasicShellDomain::StiffnessMatrix(FELinearSystem& LS, bool bsymm)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared(NE)
    for (int iel=0; iel<NE; ++iel)
    {
		FEShellElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);
        int neln = el.Nodes();
        int ndof = neln*8;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementBiphasicStiffness(el, ke, bsymm);
        
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::StiffnessMatrixSS(FELinearSystem& LS, bool bsymm)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared(NE)
    for (int iel=0; iel<NE; ++iel)
    {
		FEShellElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);
        int neln = el.Nodes();
        int ndof = neln*8;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementBiphasicStiffnessSS(el, ke, bsymm);
        
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
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
    
    const double* Mr, *Ms, *M;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    double tau = m_pMat->m_tau;
    
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
        
        // get the stress
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds c = m_pMat->Tangent(mp);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = pt.m_gradp + (pt.m_gradp - pt.m_gradpp)*(tau/dt);
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
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
                
                ke[8*i  ][8*j+6] -= kup.x;
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
                vec3d kpu = (vdotTdotv(gradMu[i], dKdE, gradMu[j])*gradp + Q*(gradMu[j]*Mu[i]))*(detJt*dt);
                vec3d kpd = (vdotTdotv(gradMu[i], dKdE, gradMd[j])*gradp + Q*(gradMd[j]*Mu[i]))*(detJt*dt);
                vec3d kqu = (vdotTdotv(gradMd[i], dKdE, gradMu[j])*gradp + Q*(gradMu[j]*Md[i]))*(detJt*dt);
                vec3d kqd = (vdotTdotv(gradMd[i], dKdE, gradMd[j])*gradp + Q*(gradMd[j]*Md[i]))*(detJt*dt);
                
                ke[8*i+6][8*j  ] -= kpu.x; ke[8*i+6][8*j+1] -= kpu.y; ke[8*i+6][8*j+2] -= kpu.z;
                
                ke[8*i+6][8*j+3] -= kpd.x; ke[8*i+6][8*j+4] -= kpd.y; ke[8*i+6][8*j+5] -= kpd.z;
                
                ke[8*i+7][8*j  ] -= kqu.x; ke[8*i+7][8*j+1] -= kqu.y; ke[8*i+7][8*j+2] -= kqu.z;
                
                ke[8*i+7][8*j+3] -= kqd.x; ke[8*i+7][8*j+4] -= kqd.y; ke[8*i+7][8*j+5] -= kqd.z;
                
            }
        
        // kpp matrix
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                double kpp = (gradMu[i]*(K*gradMu[j])*(1+tau/dt) - Phip*Mu[i]*Mu[j])*(detJt*dt);
                double kpq = (gradMu[i]*(K*gradMd[j])*(1+tau/dt) - Phip*Mu[i]*Md[j])*(detJt*dt);
                double kqp = (gradMd[i]*(K*gradMu[j])*(1+tau/dt) - Phip*Md[i]*Mu[j])*(detJt*dt);
                double kqq = (gradMd[i]*(K*gradMd[j])*(1+tau/dt) - Phip*Md[i]*Md[j])*(detJt*dt);
                
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
        
        // get the stress
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds c = m_pMat->Tangent(mp);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = pt.m_gradp;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
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
                
                ke[8*i  ][8*j+6] -= kup.x;
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
                vec3d kpu = (vdotTdotv(gradMu[i], dKdE, gradMu[j])*gradp + Q*(gradMu[j]*Mu[i]))*(detJt*dt);
                vec3d kpd = (vdotTdotv(gradMu[i], dKdE, gradMd[j])*gradp + Q*(gradMd[j]*Mu[i]))*(detJt*dt);
                vec3d kqu = (vdotTdotv(gradMd[i], dKdE, gradMu[j])*gradp + Q*(gradMu[j]*Md[i]))*(detJt*dt);
                vec3d kqd = (vdotTdotv(gradMd[i], dKdE, gradMd[j])*gradp + Q*(gradMd[j]*Md[i]))*(detJt*dt);
                
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
	FESSIShellDomain::Update(tp);

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
                if (e.DoOutput()) feLogError(e.what());
            }
        }
    }
    if (berr) throw NegativeJacobianDetected();
}

//-----------------------------------------------------------------------------
void FEBiphasicShellDomain::UpdateElementStress(int iel)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    
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
        mp.m_r0 = el.Evaluate(r0, n);
        mp.m_rt = el.Evaluate(rt, n);
        
        // get the deformation gradient and determinant
        pt.m_J = defgrad(el, pt.m_F, n);
        mat3d Fp;
        defgradp(el, Fp, n);
        mat3d Fi = pt.m_F.inverse();
        pt.m_L = (pt.m_F - Fp)*Fi / dt;

        // biphasic data
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // evaluate fluid pressure at gauss-point
        ppt.m_p = evaluate(el, pn, qn, n);
        
        // calculate the gradient of p at gauss-point
        ppt.m_gradp = gradient(el, pn, qn, n);
        
        // for biphasic materials also update the fluid flux
        //		ppt.m_w = m_pMat->Flux(mp);
        ppt.m_w = FluidFlux(mp);
        ppt.m_pa = m_pMat->Pressure(mp);
        
        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, GetFEModel()->GetTime());
        
        // calculate the solid stress at this material point
        ppt.m_ss = m_pMat->GetElasticMaterial()->Stress(mp);
        
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
        R.Assemble(el.m_node, lm, fe, true);
    }
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEBiphasicShellDomain::ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe)
{
    // get true solid and fluid densities
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
        mp.m_r0 = el.Evaluate(r0, n);
        mp.m_rt = el.Evaluate(rt, n);
        
        detJt = detJ(el, n)*gw[n];
        
        // get the force
        b = BF.force(mp);
        
        eta = el.gt(n);
        
        // evaluate apparent solid and fluid densities and mixture density
        double phiw = m_pMat->Porosity(mp);
        double rhos = (1-phiw)*m_pMat->SolidDensity(mp);
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
void FEBiphasicShellDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
    FEBiphasic* pmb = dynamic_cast<FEBiphasic*>(GetMaterial()); assert(pmb);
    
    // element stiffness matrix
    vector<int> lm;
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    for (int iel=0; iel<NE; ++iel)
    {
        FEShellElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
		FEElementMatrix ke(el);
		int neln = el.Nodes();
        int ndof = 8*neln;
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementBodyForceStiffness(bf, el, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
		ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEBiphasicShellDomain::ElementBodyForceStiffness(FEBodyForce& BF, FEShellElement &el, matrix &ke)
{
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // get true solid and fluid densities
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
        double rhos = (1-phiw)*m_pMat->SolidDensity(mp);
        double rhow = phiw*rhoTw;
        double rho = rhos + rhow;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
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
                
                kpu = (vdotTdotv(gradMu[i], dKdE, gradMu[j])*b
                       + ((b & gradMu[j]) + gradb*Mu[j])*K*gradMu[i])*(rhoTw*detJt*dt);
                kpd = (vdotTdotv(gradMu[i], dKdE, gradMd[j])*b
                       + ((b & gradMd[j]) + gradb*Md[j])*K*gradMu[i])*(rhoTw*detJt*dt);
                kqu = (vdotTdotv(gradMd[i], dKdE, gradMu[j])*b
                       + ((b & gradMu[j]) + gradb*Mu[j])*K*gradMd[i])*(rhoTw*detJt*dt);
                kqd = (vdotTdotv(gradMd[i], dKdE, gradMd[j])*b
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
    int nbf = fem.ModelLoads();
    if (nbf) {
        vec3d b(0,0,0);
        for (int i=0; i<nbf; ++i)
        {
            FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.ModelLoad(i));
            if (pbf && pbf->IsActive())
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
