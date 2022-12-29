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
#include "FEBiphasicSoluteSolidDomain.h"
#include "FECore/FEAnalysis.h"
#include "FECore/log.h"
#include <FECore/FEModel.h>
#include <FEBioMech/FEBioMech.h>
#include <FECore/FELinearSystem.h>
#include "FEBiphasicAnalysis.h"

//-----------------------------------------------------------------------------
FEBiphasicSoluteSolidDomain::FEBiphasicSoluteSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEBiphasicSoluteDomain(pfem), m_dofU(pfem), m_dofSU(pfem), m_dofR(pfem), m_dof(pfem)
{
    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
        m_dofSU.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
        m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
    }
}

//-----------------------------------------------------------------------------
//! get the material (overridden from FEDomain)
FEMaterial* FEBiphasicSoluteSolidDomain::GetMaterial()
{
	return m_pMat;
}

//-----------------------------------------------------------------------------
//! get the total dof
const FEDofList& FEBiphasicSoluteSolidDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteSolidDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
    m_pMat = dynamic_cast<FEBiphasicSolute*>(pmat);
    assert(m_pMat);
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteSolidDomain::Activate()
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
            }
        }
    }
    
    // Activate dof_P and dof_C, except when a solid element is connected to the
    // back of a shell element, in which case activate dof_Q and dof_D for those nodes.
    FEMesh& m = *GetMesh();
    for (int i=0; i<Elements(); ++i) {
        FESolidElement& el = m_Elem[i];
        int neln = el.Nodes();
        for (int j=0; j<neln; ++j)
        {
            FENode& node = m.Node(el.m_node[j]);
            if (el.m_bitfc.size()>0 && el.m_bitfc[j]) {
                node.set_active(m_dofQ);
                node.set_active(dofd  );
            }
            else {
                node.set_active(m_dofP);
                node.set_active(dofc  );
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteSolidDomain::InitMaterialPoints()
{
    FEMesh& m = *GetMesh();
    int dofc = m_dofC + m_pMat->GetSolute()->GetSoluteDOF();

    const int NE = FEElement::MAX_NODES;
    double p0[NE], c0[NE];
    
    for (int i = 0; i<(int)m_Elem.size(); ++i)
    {
        // get the solid element
        FESolidElement& el = m_Elem[i];
        
        // get the number of nodes
        int neln = el.Nodes();
        // get initial values of fluid pressure and solute concentrations
        for (int i = 0; i<neln; ++i)
        {
            FENode& node = m.Node(el.m_node[i]);
            p0[i] = node.get(m_dofP);
            c0[i] = node.get(dofc);
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
            pt.m_p = el.Evaluate(p0, n);
            pt.m_gradp = gradient(el, p0, n);
            pt.m_w = m_pMat->FluidFlux(mp);
            
            // initialize effective solute concentrations
            ps.m_c[0] = el.Evaluate(c0, n);
            ps.m_gradc[0] = gradient(el, c0, n);
            ps.m_ca[0] = m_pMat->Concentration(mp);
            ps.m_j[0] = m_pMat->SoluteFlux(mp);
            ps.m_crp[0] = pm.m_J*m_pMat->Porosity(mp)*ps.m_ca[0];
            pt.m_pa = m_pMat->Pressure(mp);
            
            // determine if solute is 'solid-bound'
            FESolute* soli = m_pMat->GetSolute();
            if (soli->m_pDiff->Diffusivity(mp).norm() == 0) ps.m_bsb[0] = true;
            
            // initialize referential solid volume fraction
            pt.m_phi0t = m_pMat->m_phi0(mp);
            
            // calculate stress
            pm.m_s = m_pMat->Stress(mp);
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data.
void FEBiphasicSoluteSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
    int dofc = m_dofC + m_pMat->GetSolute()->GetSoluteDOF();
    int dofd = m_dofD + m_pMat->GetSolute()->GetSoluteDOF();
    
    int N = el.Nodes();
    lm.resize(N*8);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = m_pMesh->Node(n);
        
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[5*i  ] = id[m_dofU[0]];
        lm[5*i+1] = id[m_dofU[1]];
        lm[5*i+2] = id[m_dofU[2]];
        
        // now the pressure dofs
        lm[5*i+3] = id[m_dofP];
        
        // concentration dofs
        lm[5*i+4] = id[dofc];
        
        // rigid rotational dofs
        // TODO: Do I really need this
        lm[5*N + 3*i  ] = id[m_dofR[0]];
        lm[5*N + 3*i+1] = id[m_dofR[1]];
        lm[5*N + 3*i+2] = id[m_dofR[2]];
    }
    
    // substitute interface dofs for solid-shell interfaces
	FESolidElement& sel = static_cast<FESolidElement&>(el);
	for (int i = 0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = m_pMesh->Node(el.m_node[i]);
            vector<int>& id = node.m_ID;
            
            // first the back-face displacement dofs
            lm[5*i  ] = id[m_dofSU[0]];
            lm[5*i+1] = id[m_dofSU[1]];
            lm[5*i+2] = id[m_dofSU[2]];
            
            // now the pressure dof (if the shell has it)
            if (id[m_dofQ] > -1) lm[5*i+3] = id[m_dofQ];
            
            // concentration dofs
            if (id[dofd] > -1) lm[5*i+4] = id[dofd];
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteSolidDomain::Reset()
{
    // reset base class
    FESolidDomain::Reset();
    
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
        ps.m_bsb.assign(nsol, false);
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
void FEBiphasicSoluteSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    int dofc = m_dofC + m_pMat->GetSolute()->GetSoluteDOF();
    int dofd = m_dofD + m_pMat->GetSolute()->GetSoluteDOF();
    
    const int NE = FEElement::MAX_NODES;
    vec3d x0[NE], xt[NE], r0, rt;
    double pn[NE], p, ct[NE], c;
    FEMesh& m = *GetMesh();
    for (size_t iel=0; iel<m_Elem.size(); ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        int neln = el.Nodes();
        for (int i=0; i<neln; ++i)
        {
            FENode& node = m.Node(el.m_node[i]);
            x0[i] = node.m_r0;
            xt[i] = node.m_rt;
            pn[i] = m.Node(el.m_node[i]).get(m_dofP);
            if (el.m_bitfc.size()>0 && el.m_bitfc[i]) {
                pn[i] = (node.m_ID[m_dofQ] != -1) ? node.get(m_dofQ) : node.get(m_dofP);
                ct[i] = (node.m_ID[dofd] != -1) ? node.get(dofd) : node.get(dofc);
            }
            else {
                pn[i] = node.get(m_dofP);
                ct[i] = node.get(dofc);
            }
        }
        
        int n = el.GaussPoints();
        for (int j=0; j<n; ++j)
        {
            r0 = el.Evaluate(x0, j);
            rt = el.Evaluate(xt, j);
            p = el.Evaluate(pn, j);
            c = el.Evaluate(ct, j);
            
            FEMaterialPoint& mp = *el.GetMaterialPoint(j);
            FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
            FEBiphasicMaterialPoint& pb = *mp.ExtractData<FEBiphasicMaterialPoint>();
            FESolutesMaterialPoint&  ps = *(mp.ExtractData<FESolutesMaterialPoint >());
            mp.m_r0 = r0;
            mp.m_rt = rt;
            
            pt.m_J = defgrad(el, pt.m_F, j);
            
            pb.m_Jp = pt.m_J;
            
            pb.m_p = p;
            pb.m_gradp = gradient(el, pn, j);
            pb.m_phi0p = pb.m_phi0t;
            ps.m_c[0] = c;
            ps.m_gradc[0] = gradient(el, ct, j);
            
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
void FEBiphasicSoluteSolidDomain::InternalForces(FEGlobalVector& R)
{
    size_t NE = m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 5*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEBiphasicSoluteSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradN;
    mat3ds s;
    
    const double* Gr, *Gs, *Gt, *H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*	gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
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
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
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
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradN = vec3d(Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i],
                          Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i],
                          Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i]);
            
            // calculate internal force
            vec3d fu = s*gradN;
            
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[5*i  ] -= fu.x*detJt;
            fe[5*i+1] -= fu.y*detJt;
            fe[5*i+2] -= fu.z*detJt;
            fe[5*i+3] -= dt*(w*gradN - divv*H[i])*detJt;
            fe[5*i+4] -= dt*(gradN*j
                             + H[i]*(crhat/J - (phiw*spt.m_ca[0] - spt.m_crp[0]/J)/dt)
                             )*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteSolidDomain::InternalForcesSS(FEGlobalVector& R)
{
    size_t NE = m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 5*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForceSS(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEBiphasicSoluteSolidDomain::ElementInternalForceSS(FESolidElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradN;
    mat3ds s;
    
    const double* Gr, *Gs, *Gt, *H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*	gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
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
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
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
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradN = vec3d(Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i],
                          Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i],
                          Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i]);
            
            // calculate internal force
            vec3d fu = s*gradN;
            
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[5*i  ] -= fu.x*detJt;
            fe[5*i+1] -= fu.y*detJt;
            fe[5*i+2] -= fu.z*detJt;
            fe[5*i+3] -= dt*(w*gradN)*detJt;
            fe[5*i+4] -= dt*(gradN*j + H[i]*(crhat/J))*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteSolidDomain::StiffnessMatrix(FELinearSystem& LS, bool bsymm)
{
    // repeat over all solid elements
    const int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FESolidElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);
        int neln = el.Nodes();
        int ndof = neln*5;
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

void FEBiphasicSoluteSolidDomain::StiffnessMatrixSS(FELinearSystem& LS, bool bsymm)
{
    // repeat over all solid elements
    const int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FESolidElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);
        int neln = el.Nodes();
        int ndof = neln*5;
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
bool FEBiphasicSoluteSolidDomain::ElementBiphasicSoluteStiffness(FESolidElement& el, matrix& ke, bool bsymm)
{
    int i, j, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(neln);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
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
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // calculate global gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
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
        vec3d vtmp,gp,gc,qpu,qcu,wc,jc;
        mat3d wu,ju;
        double qcc;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradN[i]*(s*gradN[j])) + vdotTdotv(gradN[i], C, gradN[j]))*detJ;
                ke[5*i  ][5*j  ] += Kuu[0][0]; ke[5*i  ][5*j+1] += Kuu[0][1]; ke[5*i  ][5*j+2] += Kuu[0][2];
                ke[5*i+1][5*j  ] += Kuu[1][0]; ke[5*i+1][5*j+1] += Kuu[1][1]; ke[5*i+1][5*j+2] += Kuu[1][2];
                ke[5*i+2][5*j  ] += Kuu[2][0]; ke[5*i+2][5*j+1] += Kuu[2][1]; ke[5*i+2][5*j+2] += Kuu[2][2];
                
                // calculate the kup matrix
                vtmp = -gradN[i]*H[j]*detJ;
                ke[5*i  ][5*j+3] += vtmp.x;
                ke[5*i+1][5*j+3] += vtmp.y;
                ke[5*i+2][5*j+3] += vtmp.z;
                
                // calculate the kuc matrix
                vtmp = (dTdc*gradN[i] - gradN[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*H[j]*detJ;
                ke[5*i  ][5*j+4] += vtmp.x;
                ke[5*i+1][5*j+4] += vtmp.y;
                ke[5*i+2][5*j+4] += vtmp.z;
                
                // calculate the kpu matrix
                gp = gradp+(D*gradc)*R*T*kappa/D0;
                wu = vdotTdotv(-gp, dKedE, gradN[j])
                -(((Ke*(D*gradc)) & gradN[j])*(J*dkdJ - kappa)
                  +Ke*(2*kappa*(gradN[j]*(D*gradc))))*R*T/D0
                - Ke*vdotTdotv(gradc, dDdE, gradN[j])*(kappa*R*T/D0);
                qpu = -gradN[j]*(1.0/dt);
                vtmp = (wu.transpose()*gradN[i] + qpu*H[i])*(detJ*dt);
                ke[5*i+3][5*j  ] += vtmp.x; ke[5*i+3][5*j+1] += vtmp.y; ke[5*i+3][5*j+2] += vtmp.z;
                
                // calculate the kpp matrix
                ke[5*i+3][5*j+3] -= gradN[i]*(Ke*gradN[j])*(detJ*dt);
                
                // calculate the kpc matrix
                wc = (dKedc*gp)*(-H[j])
                -Ke*((((D*(dkdc-kappa*dD0dc/D0)+dDdc*kappa)*gradc)*H[j]
                      +(D*gradN[j])*kappa)*(R*T/D0));
                ke[5*i+3][5*j+4] += (gradN[i]*wc)*(detJ*dt);
                
                // calculate the kcu matrix
                gc = -gradc*phiw + w*c/D0;
                ju = ((D*gc) & gradN[j])*(J*dkdJ)
                + vdotTdotv(gc, dDdE, gradN[j])*kappa
                + (((D*gradc) & gradN[j])*(-phis)
                   +(D*((gradN[j]*w)*2) - ((D*w) & gradN[j]))*c/D0
                   )*kappa
                +D*wu*(kappa*c/D0);
                qcu = qpu*(c*(kappa+J*phiw*dkdJ));
                vtmp = (ju.transpose()*gradN[i] + qcu*H[i])*(detJ*dt);
                ke[5*i+4][5*j  ] += vtmp.x; ke[5*i+4][5*j+1] += vtmp.y; ke[5*i+4][5*j+2] += vtmp.z;
                
                // calculate the kcp matrix
                ke[5*i+4][5*j+3] -= (gradN[i]*((D*Ke)*gradN[j]))*(kappa*c/D0)*(detJ*dt);
                
                // calculate the kcc matrix
                jc = (D*(-gradN[j]*phiw+w*(H[j]/D0)))*kappa
                +((D*dkdc+dDdc*kappa)*gc)*H[j]
                +(D*(w*(-H[j]*dD0dc/D0)+wc))*(kappa*c/D0);
                qcc = -H[j]*phiw/dt*(c*dkdc + kappa);
                ke[5*i+4][5*j+4] += (gradN[i]*jc + H[i]*qcc)*(detJ*dt);
                
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
bool FEBiphasicSoluteSolidDomain::ElementBiphasicSoluteStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm)
{
    int i, j, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(neln);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
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
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // calculate global gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds C = m_pMat->Tangent(mp);
        
        // next we get the determinant
        double J = ept.m_J;
        
        // get the fluid flux and pressure gradient
        vec3d w = ppt.m_w;
        vec3d gradp = ppt.m_gradp;
        
        // get the effective concentration, its gradient and its time derivative
        double c = spt.m_c[0];
        vec3d gradc = spt.m_gradc[0];
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->GetPermeability()->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        mat3ds dKdc = m_pMat->GetPermeability()->Tangent_Permeability_Concentration(mp, 0);
        
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
        
        // calculate all the matrices
        vec3d vtmp,gp,gc,wc,jc;
        mat3d wu,ju;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradN[i]*(s*gradN[j])) + vdotTdotv(gradN[i], C, gradN[j]))*detJ;
                
                ke[5*i  ][5*j  ] += Kuu[0][0]; ke[5*i  ][5*j+1] += Kuu[0][1]; ke[5*i  ][5*j+2] += Kuu[0][2];
                ke[5*i+1][5*j  ] += Kuu[1][0]; ke[5*i+1][5*j+1] += Kuu[1][1]; ke[5*i+1][5*j+2] += Kuu[1][2];
                ke[5*i+2][5*j  ] += Kuu[2][0]; ke[5*i+2][5*j+1] += Kuu[2][1]; ke[5*i+2][5*j+2] += Kuu[2][2];
                
                // calculate the kpu matrix
                gp = gradp+(D*gradc)*R*T*kappa/D0;
                wu = vdotTdotv(-gp, dKedE, gradN[j])
                -(((Ke*(D*gradc)) & gradN[j])*(J*dkdJ - kappa)
                  +Ke*(2*kappa*(gradN[j]*(D*gradc))))*R*T/D0
                - Ke*vdotTdotv(gradc, dDdE, gradN[j])*(kappa*R*T/D0);
                vtmp = (wu.transpose()*gradN[i])*(detJ*dt);
                ke[5*i+3][5*j  ] += vtmp.x;
                ke[5*i+3][5*j+1] += vtmp.y;
                ke[5*i+3][5*j+2] += vtmp.z;
                
                // calculate the kcu matrix
                gc = -gradc*phiw + w*c/D0;
                ju = ((D*gc) & gradN[j])*(J*dkdJ)
                + vdotTdotv(gc, dDdE, gradN[j])*kappa
                + (((D*gradc) & gradN[j])*(-phis)
                   +(D*((gradN[j]*w)*2) - ((D*w) & gradN[j]))*c/D0
                   )*kappa
                +D*wu*(kappa*c/D0);
                vtmp = (ju.transpose()*gradN[i])*(detJ*dt);
                ke[5*i+4][5*j  ] += vtmp.x;
                ke[5*i+4][5*j+1] += vtmp.y;
                ke[5*i+4][5*j+2] += vtmp.z;
                
                // calculate the kup matrix
                vtmp = -gradN[i]*H[j]*detJ;
                ke[5*i  ][5*j+3] += vtmp.x;
                ke[5*i+1][5*j+3] += vtmp.y;
                ke[5*i+2][5*j+3] += vtmp.z;
                
                // calculate the kpp matrix
                ke[5*i+3][5*j+3] -= gradN[i]*(Ke*gradN[j])*(detJ*dt);
                
                // calculate the kcp matrix
                ke[5*i+4][5*j+3] -= (gradN[i]*((D*Ke)*gradN[j]))*(kappa*c/D0)*(detJ*dt);
                
                // calculate the kuc matrix
                vtmp = (dTdc*gradN[i] - gradN[i]*(R*T*(dodc*kappa*c+osmc*dkdc*c+osmc*kappa)))*H[j]*detJ;
                ke[5*i  ][5*j+4] += vtmp.x;
                ke[5*i+1][5*j+4] += vtmp.y;
                ke[5*i+2][5*j+4] += vtmp.z;
                
                // calculate the kpc matrix
                wc = (dKedc*gp)*(-H[j])
                -Ke*((((D*(dkdc-kappa*dD0dc/D0)+dDdc*(kappa/D0))*gradc)*H[j]
                      +(D*gradN[j])*kappa)*(R*T/D0));
                ke[5*i+3][5*j+4] += (gradN[i]*wc)*(detJ*dt);
                
                // calculate the kcc matrix
                jc = (D*(-gradN[j]*phiw+w*(H[j]/D0)))*kappa
                +((D*dkdc+dDdc*kappa)*gc)*H[j]
                +(D*(w*(-H[j]*dD0dc/D0)+wc))*(kappa*c/D0);
                ke[5*i+4][5*j+4] += (gradN[i]*jc)*(detJ*dt);
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
void FEBiphasicSoluteSolidDomain::Update(const FETimeInfo& tp)
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
                if (e.DoOutput()) feLogError(e.what());
            }
        }
    }
    
    if (berr) throw NegativeJacobianDetected();
}

//-----------------------------------------------------------------------------
void FEBiphasicSoluteSolidDomain::UpdateElementStress(int iel)
{
    FEModel& fem = *GetFEModel();
    double dt = fem.GetTime().timeIncrement;
    bool sstate = (fem.GetCurrentStep()->m_nanalysis == FEBiphasicAnalysis::STEADY_STATE);
    
    int dofc = m_dofC + m_pMat->GetSolute()->GetSoluteDOF();
    int dofd = m_dofD + m_pMat->GetSolute()->GetSoluteDOF();
    
    // get the solid element
    FESolidElement& el = m_Elem[iel];
    
    // get the number of integration points
    int nint = el.GaussPoints();
    
    // get the number of nodes
    int neln = el.Nodes();
    
    // get the nodal data
    FEMesh& mesh = *m_pMesh;
    vec3d r0[FEElement::MAX_NODES];
    vec3d rt[FEElement::MAX_NODES];
    double pn[FEElement::MAX_NODES], ct[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j)
    {
        FENode& node = mesh.Node(el.m_node[j]);
        r0[j] = node.m_r0;
        rt[j] = node.m_rt;
        if (el.m_bitfc.size()>0 && el.m_bitfc[j]) {
            pn[j] = (node.m_ID[m_dofQ] != -1) ? node.get(m_dofQ) : node.get(m_dofP);
            ct[j] = (node.m_ID[dofd] != -1) ? node.get(dofd) : node.get(dofc);
        }
        else {
            pn[j] = node.get(m_dofP);
            ct[j] = node.get(dofc);
        }
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

        // solute-poroelastic data
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        
        // evaluate fluid pressure at gauss-point
        ppt.m_p = el.Evaluate(pn, n);
        
        // calculate the gradient of p at gauss-point
        ppt.m_gradp = gradient(el, pn, n);
        
        // evaluate effective solute concentration at gauss-point
        spt.m_c[0] = el.Evaluate(ct, n);
        
        // calculate the gradient of c at gauss-point
        spt.m_gradc[0] = gradient(el, ct, n);
        
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
