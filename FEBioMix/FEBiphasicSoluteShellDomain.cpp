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
#include "FEBiphasicSoluteShellDomain.h"
#include "FECore/FEMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FELinearSystem.h>
#include "FEBiphasicAnalysis.h"

//-----------------------------------------------------------------------------
FEBiphasicSoluteShellDomain::FEBiphasicSoluteShellDomain(FEModel* pfem) : FESSIShellDomain(pfem), FEBiphasicSoluteDomain(pfem), m_dof(pfem)
{
    m_dofSX = pfem->GetDOFIndex("sx");
    m_dofSY = pfem->GetDOFIndex("sy");
    m_dofSZ = pfem->GetDOFIndex("sz");
}

//-----------------------------------------------------------------------------
//! get the material (overridden from FEDomain)
FEMaterial* FEBiphasicSoluteShellDomain::GetMaterial()
{
	return m_pMat;
}

//-----------------------------------------------------------------------------
//! get the total dof
const FEDofList& FEBiphasicSoluteShellDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
    m_pMat = dynamic_cast<FEBiphasicSolute*>(pmat);
    assert(m_pMat);
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::Activate()
{
    int dofc = m_dofC + m_pMat->GetSolute()->GetSoluteDOF();
    int dofd = m_dofD + m_pMat->GetSolute()->GetSoluteDOF();
    
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
            node.set_active(dofc  );

            if (node.HasFlags(FENode::SHELL)) {
                node.set_active(m_dofQ);
                node.set_active(dofd  );
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::InitMaterialPoints()
{
    FEMesh& m = *GetMesh();
    
    const int NE = FEElement::MAX_NODES;
    double p0[NE], q0[NE], c0[NE], d0[NE];
    
    int id0 = m_pMat->GetSolute()->GetSoluteDOF();
    
    for (int i = 0; i<(int)m_Elem.size(); ++i)
    {
        // get the solid element
        FEShellElement& el = m_Elem[i];
        
        // get the number of nodes
        int neln = el.Nodes();
        // get initial values of fluid pressure and solute concentrations
        for (int i = 0; i<neln; ++i)
        {
            p0[i] = m.Node(el.m_node[i]).get(m_dofP);
            q0[i] = m.Node(el.m_node[i]).get(m_dofQ);
            c0[i] = m.Node(el.m_node[i]).get(m_dofC + id0);
            d0[i] = m.Node(el.m_node[i]).get(m_dofD + id0);
        }
        
        // get the number of integration points
        int nint = el.GaussPoints();
        
        // loop over the integration points
        for (int n = 0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEElasticMaterialPoint& pm = *(mp.ExtractData<FEElasticMaterialPoint>());
            FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            
            // initialize effective fluid pressure, its gradient, and fluid flux
            pt.m_p = evaluate(el, p0, q0, n);
            pt.m_gradp = gradient(el, p0, q0, n);
            pt.m_w = m_pMat->FluidFlux(mp);
            
            // initialize effective solute concentrations
            ps.m_c[0] = evaluate(el, c0, d0, n);
            ps.m_gradc[0] = gradient(el, c0, d0, n);
            ps.m_ca[0] = m_pMat->Concentration(mp);
            ps.m_j[0] = m_pMat->SoluteFlux(mp);
            ps.m_crp[0] = pm.m_J*m_pMat->Porosity(mp)*ps.m_ca[0];
            pt.m_pa = m_pMat->Pressure(mp);
            
            // initialize referential solid volume fraction
            pt.m_phi0t = m_pMat->m_phi0(mp);
            
            // calculate stress
            pm.m_s = m_pMat->Stress(mp);
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data.
void FEBiphasicSoluteShellDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
    int dofc = m_dofC + m_pMat->GetSolute()->GetSoluteDOF();
    int dofd = m_dofD + m_pMat->GetSolute()->GetSoluteDOF();
    int N = el.Nodes();
    int ndpn = 10;
    lm.resize(N*(ndpn+3));
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = m_pMesh->Node(n);
        
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[ndpn*i  ] = id[m_dofU[0]];
        lm[ndpn*i+1] = id[m_dofU[1]];
        lm[ndpn*i+2] = id[m_dofU[2]];
        
        // next the rotational dofs
        lm[ndpn*i+3] = id[m_dofSU[0]];
        lm[ndpn*i+4] = id[m_dofSU[1]];
        lm[ndpn*i+5] = id[m_dofSU[2]];
        
        // now the pressure dofs
        lm[ndpn*i+6] = id[m_dofP];
        lm[ndpn*i+7] = id[m_dofQ];
        
        // concentration dofs
        lm[ndpn*i+8] = id[dofc];
        lm[ndpn*i+9] = id[dofd];
        
        // rigid rotational dofs
        // TODO: Do I really need this
        lm[ndpn*N + 3*i  ] = id[m_dofR[0]];
        lm[ndpn*N + 3*i+1] = id[m_dofR[1]];
        lm[ndpn*N + 3*i+2] = id[m_dofR[2]];
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::Reset()
{
    // reset base class
    FESSIShellDomain::Reset();
    
    const int nsol = 1;
    const int nsbm = 1;

	// loop over all material points
	ForEachMaterialPoint([=](FEMaterialPoint& mp) {
        FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint&  ps = *(mp.ExtractData<FESolutesMaterialPoint >());
            
        // initialize referential solid volume fraction
        pt.m_phi0 = pt.m_phi0t = m_pMat->m_phi0(mp);
            
        // initialize multiphasic solutes
        ps.m_nsol = nsol;
        ps.m_c.assign(nsol,0);
        ps.m_ca.assign(nsol,0);
        ps.m_crp.assign(nsol, 0);
        ps.m_gradc.assign(nsol,vec3d(0,0,0));
        ps.m_k.assign(nsol, 0);
        ps.m_dkdJ.assign(nsol, 0);
        ps.m_dkdc.resize(nsol, vector<double>(nsol,0));
        ps.m_j.assign(nsol,vec3d(0,0,0));
        ps.m_nsbm = nsbm;
        ps.m_sbmr.assign(nsbm,0);
        ps.m_sbmrp.assign(nsbm,0);
        ps.m_sbmrhat.assign(nsbm,0);
        ps.m_sbmrhatp.assign(nsbm,0);
	});
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    FESSIShellDomain::PreSolveUpdate(timeInfo);
    
    int dofc = m_dofC + m_pMat->GetSolute()->GetSoluteDOF();
    int dofd = m_dofD + m_pMat->GetSolute()->GetSoluteDOF();
    
    const int NE = FEElement::MAX_NODES;
    vec3d x0[NE], xt[NE], r0, rt;
    double pn[NE], qn[NE], p;
    double cn[NE], dn[NE], c;
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
            cn[i] = m.Node(el.m_node[i]).get(dofc);
            dn[i] = m.Node(el.m_node[i]).get(dofd);
        }
        
        int n = el.GaussPoints();
        for (int j=0; j<n; ++j)
        {
            r0 = el.Evaluate(x0, j);
            rt = el.Evaluate(xt, j);
            p = evaluate(el, pn, qn, j);
            c = evaluate(el, cn, dn, j);
            
            FEMaterialPoint& mp = *el.GetMaterialPoint(j);
            FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
            FEBiphasicMaterialPoint& pb = *mp.ExtractData<FEBiphasicMaterialPoint>();
            FESolutesMaterialPoint&  ps = *(mp.ExtractData<FESolutesMaterialPoint >());
            mp.m_r0 = r0;
            mp.m_rt = rt;
            
            pt.m_J = defgrad(el, pt.m_F, j);
            
            pb.m_Jp = pt.m_J;
            
            pb.m_p = p;
            pb.m_gradp = gradient(el, pn, qn, j);
            pb.m_phi0p = pb.m_phi0t;
            
            ps.m_c[0] = c;
            ps.m_gradc[0] = gradient(el, cn, dn, j);
            
            // reset referential actual solute concentration at previous time
            ps.m_crp[0] = pt.m_J*m_pMat->Porosity(mp)*ps.m_ca[0];
            // reset referential receptor-ligand complex concentration at previous time
            ps.m_sbmrp[0] = ps.m_sbmr[0];
            // reset referential receptor-ligand complex concentration supply at previous time
            ps.m_sbmrhatp[0] = ps.m_sbmrhat[0];

            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::InternalForces(FEGlobalVector& R)
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
        int ndof = 10*el.Nodes();
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
//! calculates the internal equivalent nodal forces for solid elements

void FEBiphasicSoluteShellDomain::ElementInternalForce(FEShellElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradM, gradMu, gradMw;
    double Mu, Mw;
    mat3ds s;
    
    const double* Mr, *Ms, *M;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*	gw = el.GaussWeights();
    double eta;
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    vec3d gcnt[3];
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        
        // calculate the jacobian
        detJt = invjact(el, Ji, n);
        
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
        
        // get the solute flux
        vec3d& j = spt.m_j[0];
        
        // Evaluate porosity and solute supply and receptor-ligand kinetics
        double phiw = m_pMat->Porosity(mp);
        double crhat = 0;
        if (m_pMat->GetSolute()->m_pSupp) crhat = m_pMat->GetSolute()->m_pSupp->Supply(mp);
        
        for (i=0; i<neln; ++i)
        {
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMw = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            Mu = (1+eta)/2*M[i];
            Mw = (1-eta)/2*M[i];
            
            // calculate internal force
            vec3d fu = s*gradMu;
            vec3d fw = s*gradMw;
            
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[10*i  ] -= fu.x*detJt;
            fe[10*i+1] -= fu.y*detJt;
            fe[10*i+2] -= fu.z*detJt;
            fe[10*i+3] -= fw.x*detJt;
            fe[10*i+4] -= fw.y*detJt;
            fe[10*i+5] -= fw.z*detJt;
            fe[10*i+6] -= dt*(w*gradMu - divv*Mu)*detJt;
            fe[10*i+7] -= dt*(w*gradMw - divv*Mw)*detJt;
            fe[10*i+8] -= dt*(gradMu*j
                             + Mu*(crhat/J - (phiw*spt.m_ca[0] - spt.m_crp[0]/J)/dt)
                             )*detJt;
            fe[10*i+9] -= dt*(gradMw*j
                             + Mw*(crhat/J - (phiw*spt.m_ca[0] - spt.m_crp[0]/J)/dt)
                             )*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::InternalForcesSS(FEGlobalVector& R)
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
        int ndof = 10*el.Nodes();
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
//! calculates the internal equivalent nodal forces for solid elements

void FEBiphasicSoluteShellDomain::ElementInternalForceSS(FEShellElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradM, gradMu, gradMw;
    double Mu, Mw;
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
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        
        // calculate the jacobian
        detJt = invjact(el, Ji, n);
        
        detJt *= gw[n];
        
        // get the stress vector for this integration point
        s = pt.m_s;
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        // next we get the determinant
        double J = pt.m_J;
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        // get the solute flux
        vec3d& j = spt.m_j[0];
        
        // Evaluate solute supply and receptor-ligand kinetics
        double crhat = 0;
        if (m_pMat->GetSolute()->m_pSupp) crhat = m_pMat->GetSolute()->m_pSupp->Supply(mp);
        
        for (i=0; i<neln; ++i)
        {
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMw = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            Mu = (1+eta)/2*M[i];
            Mw = (1-eta)/2*M[i];
            
            // calculate internal force
            vec3d fu = s*gradMu;
            vec3d fw = s*gradMw;
            
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[10*i  ] -= fu.x*detJt;
            fe[10*i+1] -= fu.y*detJt;
            fe[10*i+2] -= fu.z*detJt;
            fe[10*i+3] -= fw.x*detJt;
            fe[10*i+4] -= fw.y*detJt;
            fe[10*i+5] -= fw.z*detJt;
            fe[10*i+6] -= dt*(w*gradMu)*detJt;
            fe[10*i+7] -= dt*(w*gradMw)*detJt;
            fe[10*i+8] -= dt*(gradMu*j + Mu*(crhat/J))*detJt;
            fe[10*i+9] -= dt*(gradMw*j + Mw*(crhat/J))*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::StiffnessMatrix(FELinearSystem& LS, bool bsymm)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FEShellElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);
        
        // allocate stiffness matrix
        int neln = el.Nodes();
        int ndof = neln*10;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementBiphasicSoluteStiffness(el, ke, bsymm);

		// get lm vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
}


//-----------------------------------------------------------------------------

void FEBiphasicSoluteShellDomain::StiffnessMatrixSS(FELinearSystem& LS, bool bsymm)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FEShellElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);
        int neln = el.Nodes();
        int ndof = neln*10;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementBiphasicSoluteStiffnessSS(el, ke, bsymm);

		// get lm vector
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
bool FEBiphasicSoluteShellDomain::ElementBiphasicSoluteStiffness(FEShellElement& el, matrix& ke, bool bsymm)
{
    int i, j, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    const double* Mr, *Ms, *M;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradMu(neln), gradMw(neln);
    vector<double> Mu(neln), Mw(neln);
    vec3d gradM;
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    double eta;
    
    vec3d gcnt[3];
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint&  ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint&  spt = *(mp.ExtractData<FESolutesMaterialPoint >());
        
        // calculate jacobian
        detJ = invjact(el, Ji, n)*gw[n];
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        // calculate global gradient of shape functions
        for (i=0; i<neln; ++i)
        {
            // calculate global gradient of shape functions
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu[i] = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMw[i] = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            Mu[i] = (1+eta)/2*M[i];
            Mw[i] = (1-eta)/2*M[i];
        }
        
        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds C = m_pMat->Tangent(mp);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = ppt.m_gradp;
        vec3d w = ppt.m_w;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->GetPermeability()->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        mat3ds dKdc = m_pMat->GetPermeability()->Tangent_Permeability_Concentration(mp, 0);
        
        // next we get the determinant
        double J = ept.m_J;
        
        // get the fluid flux and pressure gradient
        
        // get the effective concentration, its gradient and its time derivative
        double c = spt.m_c[0];
        vec3d gradc = spt.m_gradc[0];
        
        // evaluate the porosity and its derivative
        double phiw = m_pMat->Porosity(mp);
        double phis = 1. - phiw;
        double dpdJ = phis/J;
        
        // evaluate the solubility and its derivatives
        double kappa = spt.m_k[0];
        double dkdJ = spt.m_dkdJ[0];
        double dkdc = spt.m_dkdc[0][0];
        
        // evaluate the diffusivity tensor and its derivatives
        mat3ds D = m_pMat->GetSolute()->m_pDiff->Diffusivity(mp);
        mat3ds dDdc = m_pMat->GetSolute()->m_pDiff->Tangent_Diffusivity_Concentration(mp, 0);
        tens4dmm dDdE = m_pMat->GetSolute()->m_pDiff->Tangent_Diffusivity_Strain(mp);
        
        // evaluate the solute free diffusivity
        double D0 = m_pMat->GetSolute()->m_pDiff->Free_Diffusivity(mp);
        double dD0dc = m_pMat->GetSolute()->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,0);
        
        // evaluate the osmotic coefficient and its derivatives
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        double dodc = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp, 0);
        
        // evaluate the stress tangent with concentration
        //		mat3ds dTdc = pm->GetSolid()->Tangent_Concentration(mp, 0);
        mat3ds dTdc(0,0,0,0,0,0);
        
        // Miscellaneous constants
        mat3dd I(1);
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        
        // evaluate the effective permeability and its derivatives
        mat3ds Ki = K.inverse();
        mat3ds ImD = I-D/D0;
        mat3ds Ke = (Ki + ImD*(R*T*kappa*c/phiw/D0)).inverse();
        tens4d G = (dyad1(Ki,I) - dyad4(Ki,I)*2)*2 - ddot(dyad2(Ki,Ki),dKdE)
        +dyad1(ImD,I)*(R*T*c*J/D0/phiw*(dkdJ-kappa/phiw*dpdJ))
        +(dyad1(I,I) - dyad2(I,I)*2 - dDdE/D0)*(R*T*kappa*c/phiw/D0);
        tens4d dKedE = (dyad1(Ke,I) - 2*dyad4(Ke,I))*2 - ddot(dyad2(Ke,Ke),G);
        mat3ds Gc = -(Ki*dKdc*Ki).sym() + ImD*(R*T/phiw/D0*(dkdc*c+kappa-kappa*c/D0*dD0dc))
        +R*T*kappa*c/phiw/D0/D0*(D*dD0dc/D0 - dDdc);
        mat3ds dKedc = -(Ke*Gc*Ke).sym();
        
        // evaluate the tangents of solute supply
        double dcrhatdJ = 0;
        double dcrhatdc = 0;
        if (m_pMat->GetSolute()->m_pSupp)
        {
            dcrhatdJ = m_pMat->GetSolute()->m_pSupp->Tangent_Supply_Strain(mp);
            double dcrhatdcr = m_pMat->GetSolute()->m_pSupp->Tangent_Supply_Concentration(mp);
            dcrhatdc = J*phiw*(kappa + c*dkdc)*dcrhatdcr;
        }
        
        // calculate all the matrices
        vec3d vtmp,gp,gc,qpu,qpw,qcu,qcw,wc,wd,jc,jd;
        mat3d wu,ww,ju,jw;
        double qcc, qcd;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradMu[i]*(s*gradMu[j])) + vdotTdotv(gradMu[i], C, gradMu[j]))*detJ;
                mat3d Kuw = (mat3dd(gradMu[i]*(s*gradMw[j])) + vdotTdotv(gradMu[i], C, gradMw[j]))*detJ;
                mat3d Kwu = (mat3dd(gradMw[i]*(s*gradMu[j])) + vdotTdotv(gradMw[i], C, gradMu[j]))*detJ;
                mat3d Kww = (mat3dd(gradMw[i]*(s*gradMw[j])) + vdotTdotv(gradMw[i], C, gradMw[j]))*detJ;
                
                ke[10*i  ][10*j  ] += Kuu[0][0]; ke[10*i  ][10*j+1] += Kuu[0][1]; ke[10*i  ][10*j+2] += Kuu[0][2];
                ke[10*i+1][10*j  ] += Kuu[1][0]; ke[10*i+1][10*j+1] += Kuu[1][1]; ke[10*i+1][10*j+2] += Kuu[1][2];
                ke[10*i+2][10*j  ] += Kuu[2][0]; ke[10*i+2][10*j+1] += Kuu[2][1]; ke[10*i+2][10*j+2] += Kuu[2][2];
                
                ke[10*i  ][10*j+3] += Kuw[0][0]; ke[10*i  ][10*j+4] += Kuw[0][1]; ke[10*i  ][10*j+5] += Kuw[0][2];
                ke[10*i+1][10*j+3] += Kuw[1][0]; ke[10*i+1][10*j+4] += Kuw[1][1]; ke[10*i+1][10*j+5] += Kuw[1][2];
                ke[10*i+2][10*j+3] += Kuw[2][0]; ke[10*i+2][10*j+4] += Kuw[2][1]; ke[10*i+2][10*j+5] += Kuw[2][2];
                
                ke[10*i+3][10*j  ] += Kwu[0][0]; ke[10*i+3][10*j+1] += Kwu[0][1]; ke[10*i+3][10*j+2] += Kwu[0][2];
                ke[10*i+4][10*j  ] += Kwu[1][0]; ke[10*i+4][10*j+1] += Kwu[1][1]; ke[10*i+4][10*j+2] += Kwu[1][2];
                ke[10*i+5][10*j  ] += Kwu[2][0]; ke[10*i+5][10*j+1] += Kwu[2][1]; ke[10*i+5][10*j+2] += Kwu[2][2];
                
                ke[10*i+3][10*j+3] += Kww[0][0]; ke[10*i+3][10*j+4] += Kww[0][1]; ke[10*i+3][10*j+5] += Kww[0][2];
                ke[10*i+4][10*j+3] += Kww[1][0]; ke[10*i+4][10*j+4] += Kww[1][1]; ke[10*i+4][10*j+5] += Kww[1][2];
                ke[10*i+5][10*j+3] += Kww[2][0]; ke[10*i+5][10*j+4] += Kww[2][1]; ke[10*i+5][10*j+5] += Kww[2][2];
                
                // calculate the kup matrix
                vec3d kup = gradMu[i]*(-Mu[j]*detJ);
                vec3d kuq = gradMu[i]*(-Mw[j]*detJ);
                vec3d kwp = gradMw[i]*(-Mu[j]*detJ);
                vec3d kwq = gradMw[i]*(-Mw[j]*detJ);
                
                ke[10*i  ][10*j+6] += kup.x; ke[10*i  ][10*j+7] += kuq.x;
                ke[10*i+1][10*j+6] += kup.y; ke[10*i+1][10*j+7] += kuq.y;
                ke[10*i+2][10*j+6] += kup.z; ke[10*i+2][10*j+7] += kuq.z;
                
                ke[10*i+3][10*j+6] += kwp.x; ke[10*i+3][10*j+7] += kwq.x;
                ke[10*i+4][10*j+6] += kwp.y; ke[10*i+4][10*j+7] += kwq.y;
                ke[10*i+5][10*j+6] += kwp.z; ke[10*i+5][10*j+7] += kwq.z;
                
                // calculate the kuc matrix
                vec3d kuc = (dTdc*gradMu[i] - gradMu[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*Mu[j]*detJ;
                vec3d kud = (dTdc*gradMu[i] - gradMu[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*Mw[j]*detJ;
                vec3d kwc = (dTdc*gradMw[i] - gradMw[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*Mu[j]*detJ;
                vec3d kwd = (dTdc*gradMw[i] - gradMw[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*Mw[j]*detJ;
                
                ke[10*i  ][10*j+8] += kuc.x; ke[10*i  ][10*j+9] += kud.x;
                ke[10*i+1][10*j+8] += kuc.y; ke[10*i+1][10*j+9] += kud.y;
                ke[10*i+2][10*j+8] += kuc.z; ke[10*i+2][10*j+9] += kud.z;
                
                ke[10*i+3][10*j+8] += kwc.x; ke[10*i+3][10*j+9] += kwd.x;
                ke[10*i+4][10*j+8] += kwc.y; ke[10*i+4][10*j+9] += kwd.y;
                ke[10*i+5][10*j+8] += kwc.z; ke[10*i+5][10*j+9] += kwd.z;
                
                // calculate the kpu matrix
                gp = gradp+(D*gradc)*R*T*kappa/D0;
                wu = vdotTdotv(-gp, dKedE, gradMu[j])
                -(((Ke*(D*gradc)) & gradMu[j])*(J*dkdJ - kappa)
                  +Ke*(2*kappa*(gradMu[j]*(D*gradc))))*R*T/D0
                - Ke*vdotTdotv(gradc, dDdE, gradMu[j])*(kappa*R*T/D0);
                ww = vdotTdotv(-gp, dKedE, gradMw[j])
                -(((Ke*(D*gradc)) & gradMw[j])*(J*dkdJ - kappa)
                  +Ke*(2*kappa*(gradMw[j]*(D*gradc))))*R*T/D0
                - Ke*vdotTdotv(gradc, dDdE, gradMw[j])*(kappa*R*T/D0);
                qpu = -gradMu[j]*(1.0/dt);
                qpw = -gradMw[j]*(1.0/dt);
                vec3d kpu = (wu.transpose()*gradMu[i] + qpu*Mu[i])*(detJ*dt);
                vec3d kpw = (ww.transpose()*gradMu[i] + qpw*Mu[i])*(detJ*dt);
                vec3d kqu = (wu.transpose()*gradMw[i] + qpu*Mw[i])*(detJ*dt);
                vec3d kqw = (ww.transpose()*gradMw[i] + qpw*Mw[i])*(detJ*dt);
                ke[10*i+6][10*j  ] += kpu.x; ke[10*i+6][10*j+1] += kpu.y; ke[10*i+6][10*j+2] += kpu.z;
                ke[10*i+6][10*j+3] += kpw.x; ke[10*i+6][10*j+4] += kpw.y; ke[10*i+6][10*j+5] += kpw.z;
                ke[10*i+7][10*j  ] += kqu.x; ke[10*i+7][10*j+1] += kqu.y; ke[10*i+7][10*j+2] += kqu.z;
                ke[10*i+7][10*j+3] += kqw.x; ke[10*i+7][10*j+4] += kqw.y; ke[10*i+7][10*j+5] += kqw.z;
                
                // calculate the kpp matrix
                ke[10*i+6][10*j+6] -= gradMu[i]*(Ke*gradMu[j])*(detJ*dt); ke[10*i+6][10*j+7] -= gradMu[i]*(Ke*gradMw[j])*(detJ*dt);
                ke[10*i+7][10*j+6] -= gradMw[i]*(Ke*gradMu[j])*(detJ*dt); ke[10*i+7][10*j+7] -= gradMw[i]*(Ke*gradMw[j])*(detJ*dt);
                
                // calculate the kpc matrix
                wc = (dKedc*gp)*(-Mu[j])
                -Ke*((((D*(dkdc-kappa*dD0dc/D0)+dDdc*kappa)*gradc)*Mu[j]
                      +(D*gradMu[j])*kappa)*(R*T/D0));
                wd = (dKedc*gp)*(-Mw[j])
                -Ke*((((D*(dkdc-kappa*dD0dc/D0)+dDdc*kappa)*gradc)*Mw[j]
                      +(D*gradMw[j])*kappa)*(R*T/D0));
                ke[10*i+6][10*j+8] += (gradMu[i]*wc)*(detJ*dt);
                ke[10*i+6][10*j+9] += (gradMu[i]*wd)*(detJ*dt);
                ke[10*i+7][10*j+8] += (gradMw[i]*wc)*(detJ*dt);
                ke[10*i+7][10*j+9] += (gradMw[i]*wd)*(detJ*dt);
                
                // calculate the kcu matrix
                gc = -gradc*phiw + w*c/D0;
                ju = ((D*gc) & gradMu[j])*(J*dkdJ)
                + vdotTdotv(gc, dDdE, gradMu[j])*kappa
                + (((D*gradc) & gradMu[j])*(-phis)
                   +(D*((gradMu[j]*w)*2) - ((D*w) & gradMu[j]))*c/D0
                   )*kappa
                +D*wu*(kappa*c/D0);
                jw = ((D*gc) & gradMw[j])*(J*dkdJ)
                + vdotTdotv(gc, dDdE, gradMw[j])*kappa
                + (((D*gradc) & gradMw[j])*(-phis)
                   +(D*((gradMw[j]*w)*2) - ((D*w) & gradMw[j]))*c/D0
                   )*kappa
                +D*ww*(kappa*c/D0);
                qcu = qpu*(c*(kappa+J*phiw*dkdJ));
                qcw = qpw*(c*(kappa+J*phiw*dkdJ));
                vec3d kcu = (ju.transpose()*gradMu[i] + qcu*Mu[i])*(detJ*dt);
                vec3d kcw = (jw.transpose()*gradMu[i] + qcw*Mu[i])*(detJ*dt);
                vec3d kdu = (ju.transpose()*gradMw[i] + qcu*Mw[i])*(detJ*dt);
                vec3d kdw = (jw.transpose()*gradMw[i] + qcw*Mw[i])*(detJ*dt);
                ke[10*i+8][10*j  ] += kcu.x; ke[10*i+8][10*j+1] += kcu.y; ke[10*i+8][10*j+2] += kcu.z;
                ke[10*i+8][10*j+3] += kcw.x; ke[10*i+8][10*j+4] += kcw.y; ke[10*i+8][10*j+5] += kcw.z;
                ke[10*i+9][10*j  ] += kdu.x; ke[10*i+9][10*j+1] += kdu.y; ke[10*i+9][10*j+2] += kdu.z;
                ke[10*i+9][10*j+3] += kdw.x; ke[10*i+9][10*j+4] += kdw.y; ke[10*i+9][10*j+5] += kdw.z;
                
                // calculate the kcp matrix
                ke[10*i+8][10*j+6] -= (gradMu[i]*((D*Ke)*gradMu[j]))*(kappa*c/D0)*(detJ*dt); ke[10*i+8][10*j+7] -= (gradMu[i]*((D*Ke)*gradMw[j]))*(kappa*c/D0)*(detJ*dt);
                ke[10*i+9][10*j+6] -= (gradMw[i]*((D*Ke)*gradMu[j]))*(kappa*c/D0)*(detJ*dt); ke[10*i+9][10*j+7] -= (gradMw[i]*((D*Ke)*gradMw[j]))*(kappa*c/D0)*(detJ*dt);
                
                // calculate the kcc matrix
                jc = (D*(-gradMu[j]*phiw+w*(Mu[j]/D0)))*kappa
                +((D*dkdc+dDdc*kappa)*gc)*Mu[j]
                +(D*(w*(-Mu[j]*dD0dc/D0)+wc))*(kappa*c/D0);
                jd = (D*(-gradMw[j]*phiw+w*(Mw[j]/D0)))*kappa
                +((D*dkdc+dDdc*kappa)*gc)*Mw[j]
                +(D*(w*(-Mw[j]*dD0dc/D0)+wd))*(kappa*c/D0);
                qcc = -Mu[j]*phiw/dt*(c*dkdc + kappa);
                qcd = -Mw[j]*phiw/dt*(c*dkdc + kappa);
                ke[10*i+8][10*j+8] += (gradMu[i]*jc + Mu[i]*qcc)*(detJ*dt);
                ke[10*i+8][10*j+9] += (gradMu[i]*jd + Mu[i]*qcd)*(detJ*dt);
                ke[10*i+9][10*j+8] += (gradMw[i]*jc + Mw[i]*qcc)*(detJ*dt);
                ke[10*i+9][10*j+9] += (gradMw[i]*jd + Mw[i]*qcd)*(detJ*dt);
                
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
    if (bsymm) {
        for (i=0; i<5*neln; ++i)
            for (j=i+1; j<5*neln; ++j) {
                tmp = 0.5*(ke[i][j]+ke[j][i]);
                ke[i][j] = ke[j][i] = tmp;
            }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//! for steady-state response (zero solid velocity, zero time derivative of
//! solute concentration)
//!
bool FEBiphasicSoluteShellDomain::ElementBiphasicSoluteStiffnessSS(FEShellElement& el, matrix& ke, bool bsymm)
{
    int i, j, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    const double* Mr, *Ms, *M;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradMu(neln), gradMw(neln);
    vector<double> Mu(neln), Mw(neln);
    vec3d gradM;
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    double eta;
    
    vec3d gcnt[3];
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint&  ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint&  spt = *(mp.ExtractData<FESolutesMaterialPoint >());
        
        // calculate jacobian
        detJ = invjact(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        // calculate global gradient of shape functions
        for (i=0; i<neln; ++i)
        {
            // calculate global gradient of shape functions
            gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu[i] = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMw[i] = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            Mu[i] = (1+eta)/2*M[i];
            Mw[i] = (1-eta)/2*M[i];
        }
        
        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds C = m_pMat->Tangent(mp);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = ppt.m_gradp;
        vec3d w = ppt.m_w;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->GetPermeability()->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        mat3ds dKdc = m_pMat->GetPermeability()->Tangent_Permeability_Concentration(mp, 0);
        
        // next we get the determinant
        double J = ept.m_J;
        
        // get the fluid flux and pressure gradient
        
        // get the effective concentration, its gradient and its time derivative
        double c = spt.m_c[0];
        vec3d gradc = spt.m_gradc[0];
        
        // evaluate the porosity and its derivative
        double phiw = m_pMat->Porosity(mp);
        double phis = 1. - phiw;
        double dpdJ = phis/J;
        
        // evaluate the solubility and its derivatives
        double kappa = spt.m_k[0];
        double dkdJ = spt.m_dkdJ[0];
        double dkdc = spt.m_dkdc[0][0];
        
        // evaluate the diffusivity tensor and its derivatives
        mat3ds D = m_pMat->GetSolute()->m_pDiff->Diffusivity(mp);
        mat3ds dDdc = m_pMat->GetSolute()->m_pDiff->Tangent_Diffusivity_Concentration(mp, 0);
        tens4dmm dDdE = m_pMat->GetSolute()->m_pDiff->Tangent_Diffusivity_Strain(mp);
        
        // evaluate the solute free diffusivity
        double D0 = m_pMat->GetSolute()->m_pDiff->Free_Diffusivity(mp);
        double dD0dc = m_pMat->GetSolute()->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,0);
        
        // evaluate the osmotic coefficient and its derivatives
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        double dodc = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp, 0);
        
        // evaluate the stress tangent with concentration
        //		mat3ds dTdc = pm->GetSolid()->Tangent_Concentration(mp, 0);
        mat3ds dTdc(0,0,0,0,0,0);
        
        // Miscellaneous constants
        mat3dd I(1);
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        
        // evaluate the effective permeability and its derivatives
        mat3ds Ki = K.inverse();
        mat3ds ImD = I-D/D0;
        mat3ds Ke = (Ki + ImD*(R*T*kappa*c/phiw/D0)).inverse();
        tens4d G = (dyad1(Ki,I) - dyad4(Ki,I)*2)*2 - ddot(dyad2(Ki,Ki),dKdE)
        +dyad1(ImD,I)*(R*T*c*J/D0/phiw*(dkdJ-kappa/phiw*dpdJ))
        +(dyad1(I,I) - dyad2(I,I)*2 - dDdE/D0)*(R*T*kappa*c/phiw/D0);
        tens4d dKedE = (dyad1(Ke,I) - 2*dyad4(Ke,I))*2 - ddot(dyad2(Ke,Ke),G);
        mat3ds Gc = -(Ki*dKdc*Ki).sym() + ImD*(R*T/phiw/D0*(dkdc*c+kappa-kappa*c/D0*dD0dc))
        +R*T*kappa*c/phiw/D0/D0*(D*dD0dc/D0 - dDdc);
        mat3ds dKedc = -(Ke*Gc*Ke).sym();
        
        // evaluate the tangents of solute supply
        double dcrhatdJ = 0;
        double dcrhatdc = 0;
        if (m_pMat->GetSolute()->m_pSupp)
        {
            dcrhatdJ = m_pMat->GetSolute()->m_pSupp->Tangent_Supply_Strain(mp);
            double dcrhatdcr = m_pMat->GetSolute()->m_pSupp->Tangent_Supply_Concentration(mp);
            dcrhatdc = J*phiw*(kappa + c*dkdc)*dcrhatdcr;
        }
        
        // calculate all the matrices
        vec3d vtmp,gp,gc,wc,wd,jc,jd;
        mat3d wu,ww,ju,jw;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradMu[i]*(s*gradMu[j])) + vdotTdotv(gradMu[i], C, gradMu[j]))*detJ;
                mat3d Kuw = (mat3dd(gradMu[i]*(s*gradMw[j])) + vdotTdotv(gradMu[i], C, gradMw[j]))*detJ;
                mat3d Kwu = (mat3dd(gradMw[i]*(s*gradMu[j])) + vdotTdotv(gradMw[i], C, gradMu[j]))*detJ;
                mat3d Kww = (mat3dd(gradMw[i]*(s*gradMw[j])) + vdotTdotv(gradMw[i], C, gradMw[j]))*detJ;
                
                ke[10*i  ][10*j  ] += Kuu[0][0]; ke[10*i  ][10*j+1] += Kuu[0][1]; ke[10*i  ][10*j+2] += Kuu[0][2];
                ke[10*i+1][10*j  ] += Kuu[1][0]; ke[10*i+1][10*j+1] += Kuu[1][1]; ke[10*i+1][10*j+2] += Kuu[1][2];
                ke[10*i+2][10*j  ] += Kuu[2][0]; ke[10*i+2][10*j+1] += Kuu[2][1]; ke[10*i+2][10*j+2] += Kuu[2][2];
                
                ke[10*i  ][10*j+3] += Kuw[0][0]; ke[10*i  ][10*j+4] += Kuw[0][1]; ke[10*i  ][10*j+5] += Kuw[0][2];
                ke[10*i+1][10*j+3] += Kuw[1][0]; ke[10*i+1][10*j+4] += Kuw[1][1]; ke[10*i+1][10*j+5] += Kuw[1][2];
                ke[10*i+2][10*j+3] += Kuw[2][0]; ke[10*i+2][10*j+4] += Kuw[2][1]; ke[10*i+2][10*j+5] += Kuw[2][2];
                
                ke[10*i+3][10*j  ] += Kwu[0][0]; ke[10*i+3][10*j+1] += Kwu[0][1]; ke[10*i+3][10*j+2] += Kwu[0][2];
                ke[10*i+4][10*j  ] += Kwu[1][0]; ke[10*i+4][10*j+1] += Kwu[1][1]; ke[10*i+4][10*j+2] += Kwu[1][2];
                ke[10*i+5][10*j  ] += Kwu[2][0]; ke[10*i+5][10*j+1] += Kwu[2][1]; ke[10*i+5][10*j+2] += Kwu[2][2];
                
                ke[10*i+3][10*j+3] += Kww[0][0]; ke[10*i+3][10*j+4] += Kww[0][1]; ke[10*i+3][10*j+5] += Kww[0][2];
                ke[10*i+4][10*j+3] += Kww[1][0]; ke[10*i+4][10*j+4] += Kww[1][1]; ke[10*i+4][10*j+5] += Kww[1][2];
                ke[10*i+5][10*j+3] += Kww[2][0]; ke[10*i+5][10*j+4] += Kww[2][1]; ke[10*i+5][10*j+5] += Kww[2][2];
                
                // calculate the kup matrix
                vec3d kup = gradMu[i]*(-Mu[j]*detJ);
                vec3d kuq = gradMu[i]*(-Mw[j]*detJ);
                vec3d kwp = gradMw[i]*(-Mu[j]*detJ);
                vec3d kwq = gradMw[i]*(-Mw[j]*detJ);
                
                ke[10*i  ][10*j+6] += kup.z; ke[10*i  ][10*j+7] += kuq.x;
                ke[10*i+1][10*j+6] += kup.y; ke[10*i+1][10*j+7] += kuq.y;
                ke[10*i+2][10*j+6] += kup.z; ke[10*i+2][10*j+7] += kuq.z;
                
                ke[10*i+3][10*j+6] += kwp.x; ke[10*i+3][10*j+7] += kwq.x;
                ke[10*i+4][10*j+6] += kwp.y; ke[10*i+4][10*j+7] += kwq.y;
                ke[10*i+5][10*j+6] += kwp.z; ke[10*i+5][10*j+7] += kwq.z;
                
                // calculate the kuc matrix
                vec3d kuc = (dTdc*gradMu[i] - gradMu[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*Mu[j]*detJ;
                vec3d kud = (dTdc*gradMu[i] - gradMu[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*Mw[j]*detJ;
                vec3d kwc = (dTdc*gradMw[i] - gradMw[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*Mu[j]*detJ;
                vec3d kwd = (dTdc*gradMw[i] - gradMw[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*Mw[j]*detJ;
                
                ke[10*i  ][10*j+8] += kuc.z; ke[10*i  ][10*j+9] += kud.x;
                ke[10*i+1][10*j+8] += kuc.y; ke[10*i+1][10*j+9] += kud.y;
                ke[10*i+2][10*j+8] += kuc.z; ke[10*i+2][10*j+9] += kud.z;
                
                ke[10*i+3][10*j+8] += kwc.x; ke[10*i+3][10*j+9] += kwd.x;
                ke[10*i+4][10*j+8] += kwc.y; ke[10*i+4][10*j+9] += kwd.y;
                ke[10*i+5][10*j+8] += kwc.z; ke[10*i+5][10*j+9] += kwd.z;
                
                // calculate the kpu matrix
                gp = gradp+(D*gradc)*R*T*kappa/D0;
                wu = vdotTdotv(-gp, dKedE, gradMu[j])
                -(((Ke*(D*gradc)) & gradMu[j])*(J*dkdJ - kappa)
                  +Ke*(2*kappa*(gradMu[j]*(D*gradc))))*R*T/D0
                - Ke*vdotTdotv(gradc, dDdE, gradMu[j])*(kappa*R*T/D0);
                ww = vdotTdotv(-gp, dKedE, gradMw[j])
                -(((Ke*(D*gradc)) & gradMw[j])*(J*dkdJ - kappa)
                  +Ke*(2*kappa*(gradMw[j]*(D*gradc))))*R*T/D0
                - Ke*vdotTdotv(gradc, dDdE, gradMw[j])*(kappa*R*T/D0);
                vec3d kpu = (wu.transpose()*gradMu[i])*(detJ*dt);
                vec3d kpw = (ww.transpose()*gradMu[i])*(detJ*dt);
                vec3d kqu = (wu.transpose()*gradMw[i])*(detJ*dt);
                vec3d kqw = (ww.transpose()*gradMw[i])*(detJ*dt);
                ke[10*i+6][10*j  ] += kpu.x; ke[10*i+6][10*j+1] += kpu.y; ke[10*i+6][10*j+2] += kpu.z;
                ke[10*i+6][10*j+3] += kpw.x; ke[10*i+6][10*j+4] += kpw.y; ke[10*i+6][10*j+5] += kpw.z;
                ke[10*i+7][10*j  ] += kqu.x; ke[10*i+7][10*j+1] += kqu.y; ke[10*i+7][10*j+2] += kqu.z;
                ke[10*i+7][10*j+3] += kqw.x; ke[10*i+7][10*j+4] += kqw.y; ke[10*i+7][10*j+5] += kqw.z;
                
                // calculate the kpp matrix
                ke[10*i+6][10*j+6] -= gradMu[i]*(Ke*gradMu[j])*(detJ*dt); ke[10*i+6][10*j+7] -= gradMu[i]*(Ke*gradMw[j])*(detJ*dt);
                ke[10*i+7][10*j+6] -= gradMw[i]*(Ke*gradMu[j])*(detJ*dt); ke[10*i+7][10*j+7] -= gradMw[i]*(Ke*gradMw[j])*(detJ*dt);
                
                // calculate the kpc matrix
                wc = (dKedc*gp)*(-Mu[j])
                -Ke*((((D*(dkdc-kappa*dD0dc/D0)+dDdc*kappa)*gradc)*Mu[j]
                      +(D*gradMu[j])*kappa)*(R*T/D0));
                wd = (dKedc*gp)*(-Mw[j])
                -Ke*((((D*(dkdc-kappa*dD0dc/D0)+dDdc*kappa)*gradc)*Mw[j]
                      +(D*gradMw[j])*kappa)*(R*T/D0));
                ke[10*i+6][10*j+8] += (gradMu[i]*wc)*(detJ*dt);
                ke[10*i+6][10*j+9] += (gradMu[i]*wd)*(detJ*dt);
                ke[10*i+7][10*j+8] += (gradMw[i]*wc)*(detJ*dt);
                ke[10*i+7][10*j+9] += (gradMw[i]*wd)*(detJ*dt);
                
                // calculate the kcu matrix
                gc = -gradc*phiw + w*c/D0;
                ju = ((D*gc) & gradMu[j])*(J*dkdJ)
                + vdotTdotv(gc, dDdE, gradMu[j])*kappa
                + (((D*gradc) & gradMu[j])*(-phis)
                   +(D*((gradMu[j]*w)*2) - ((D*w) & gradMu[j]))*c/D0
                   )*kappa
                +D*wu*(kappa*c/D0);
                jw = ((D*gc) & gradMw[j])*(J*dkdJ)
                + vdotTdotv(gc, dDdE, gradMw[j])*kappa
                + (((D*gradc) & gradMw[j])*(-phis)
                   +(D*((gradMw[j]*w)*2) - ((D*w) & gradMw[j]))*c/D0
                   )*kappa
                +D*ww*(kappa*c/D0);
                vec3d kcu = (ju.transpose()*gradMu[i])*(detJ*dt);
                vec3d kcw = (jw.transpose()*gradMu[i])*(detJ*dt);
                vec3d kdu = (ju.transpose()*gradMw[i])*(detJ*dt);
                vec3d kdw = (jw.transpose()*gradMw[i])*(detJ*dt);
                ke[10*i+8][10*j  ] += kcu.x; ke[10*i+8][10*j+1] += kcu.y; ke[10*i+8][10*j+2] += kcu.z;
                ke[10*i+8][10*j+3] += kcw.x; ke[10*i+8][10*j+4] += kcw.y; ke[10*i+8][10*j+5] += kcw.z;
                ke[10*i+9][10*j  ] += kdu.x; ke[10*i+9][10*j+1] += kdu.y; ke[10*i+9][10*j+2] += kdu.z;
                ke[10*i+9][10*j+3] += kdw.x; ke[10*i+9][10*j+4] += kdw.y; ke[10*i+9][10*j+5] += kdw.z;
                
                // calculate the kcp matrix
                ke[10*i+8][10*j+6] -= (gradMu[i]*((D*Ke)*gradMu[j]))*(kappa*c/D0)*(detJ*dt); ke[10*i+8][10*j+7] -= (gradMu[i]*((D*Ke)*gradMw[j]))*(kappa*c/D0)*(detJ*dt);
                ke[10*i+9][10*j+6] -= (gradMw[i]*((D*Ke)*gradMu[j]))*(kappa*c/D0)*(detJ*dt); ke[10*i+9][10*j+7] -= (gradMw[i]*((D*Ke)*gradMw[j]))*(kappa*c/D0)*(detJ*dt);
                
                // calculate the kcc matrix
                jc = (D*(-gradMu[j]*phiw+w*(Mu[j]/D0)))*kappa
                +((D*dkdc+dDdc*kappa)*gc)*Mu[j]
                +(D*(w*(-Mu[j]*dD0dc/D0)+wc))*(kappa*c/D0);
                jd = (D*(-gradMw[j]*phiw+w*(Mw[j]/D0)))*kappa
                +((D*dkdc+dDdc*kappa)*gc)*Mw[j]
                +(D*(w*(-Mw[j]*dD0dc/D0)+wd))*(kappa*c/D0);
                ke[10*i+8][10*j+8] += (gradMu[i]*jc)*(detJ*dt);
                ke[10*i+8][10*j+9] += (gradMu[i]*jd)*(detJ*dt);
                ke[10*i+9][10*j+8] += (gradMw[i]*jc)*(detJ*dt);
                ke[10*i+9][10*j+9] += (gradMw[i]*jd)*(detJ*dt);
                
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
    if (bsymm) {
        for (i=0; i<5*neln; ++i)
            for (j=i+1; j<5*neln; ++j) {
                tmp = 0.5*(ke[i][j]+ke[j][i]);
                ke[i][j] = ke[j][i] = tmp;
            }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteShellDomain::Update(const FETimeInfo& tp)
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
void FEBiphasicSoluteShellDomain::UpdateElementStress(int iel)
{
    FEModel& fem = *GetFEModel();
    double dt = fem.GetTime().timeIncrement;
    bool sstate = (fem.GetCurrentStep()->m_nanalysis == FEBiphasicAnalysis::STEADY_STATE);
    
    // get the solid element
    FEShellElement& el = m_Elem[iel];
    
    // get the number of integration points
    int nint = el.GaussPoints();
    
    // get the number of nodes
    int neln = el.Nodes();
    
    // get the biphasic-solute material
    int id0 = m_pMat->GetSolute()->GetSoluteDOF();
    
    // get the nodal data
    FEMesh& mesh = *m_pMesh;
    vec3d r0[FEElement::MAX_NODES];
    vec3d rt[FEElement::MAX_NODES];
    double pn[FEElement::MAX_NODES], qn[FEElement::MAX_NODES];
    double cn[FEElement::MAX_NODES], dn[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j)
    {
        r0[j] = mesh.Node(el.m_node[j]).m_r0;
        rt[j] = mesh.Node(el.m_node[j]).m_rt;
        pn[j] = mesh.Node(el.m_node[j]).get(m_dofP);
        qn[j] = mesh.Node(el.m_node[j]).get(m_dofQ);
        cn[j] = mesh.Node(el.m_node[j]).get(m_dofC + id0);
        dn[j] = mesh.Node(el.m_node[j]).get(m_dofD + id0);
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

        // biphasic-solute data
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        
        // evaluate fluid pressure at gauss-point
        ppt.m_p = evaluate(el, pn, qn, n);
        
        // calculate the gradient of p at gauss-point
        ppt.m_gradp = gradient(el, pn, qn, n);
        
        // evaluate effective solute concentration at gauss-point
        spt.m_c[0] = evaluate(el, cn, dn, n);
        
        // calculate the gradient of c at gauss-point
        spt.m_gradc[0] = gradient(el, cn, dn, n);
        
        // for biphasic-solute materials also update the porosity, fluid and solute fluxes
        // and evaluate the actual fluid pressure and solute concentration
        ppt.m_w = m_pMat->FluidFlux(mp);
        ppt.m_pa = m_pMat->Pressure(mp);
        spt.m_j[0] = m_pMat->SoluteFlux(mp);
        spt.m_ca[0] = m_pMat->Concentration(mp);
        if (m_pMat->GetSolute()->m_pSupp)
        {
            if (sstate)
                spt.m_sbmr[0] = m_pMat->GetSolute()->m_pSupp->ReceptorLigandConcentrationSS(mp);
            else {
                // update m_crc using midpoint rule
                spt.m_sbmrhat[0] = m_pMat->GetSolute()->m_pSupp->ReceptorLigandSupply(mp);
                spt.m_sbmr[0] = spt.m_sbmrp[0] + (spt.m_sbmrhat[0]+spt.m_sbmrhatp[0])/2*dt;
                // update phi0 using backward difference integration
                
                // NOTE: MolarMass was removed since not used
                ppt.m_phi0hat = 0;
                //				ppt.m_phi0hat = pmb->GetSolid()->MolarMass()/pmb->GetSolid()->Density()*pmb->GetSolute()->m_pSupp->SolidSupply(mp);
                
                ppt.m_phi0t = ppt.m_phi0p + ppt.m_phi0hat*dt;
            }
        }
        
        m_pMat->PartitionCoefficientFunctions(mp, spt.m_k[0], spt.m_dkdJ[0], spt.m_dkdc[0][0]);

        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, GetFEModel()->GetTime());
        
        // calculate the solid stress at this material point
        ppt.m_ss = m_pMat->GetElasticMaterial()->Stress(mp);
        
        // calculate the stress at this material point (must be done after evaluating m_pa)
        pt.m_s = m_pMat->Stress(mp);
    }
}
