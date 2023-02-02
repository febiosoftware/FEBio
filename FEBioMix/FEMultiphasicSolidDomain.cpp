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
#include "FEMultiphasicSolidDomain.h"
#include "FEMultiphasicMultigeneration.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FECore/DOFS.h>
#include <FEBioMech/FEBioMech.h>
#include <FECore/FELinearSystem.h>
#include <FECore/sys.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
FEMultiphasicSolidDomain::FEMultiphasicSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEMultiphasicDomain(pfem), m_dofU(pfem), m_dofSU(pfem), m_dofR(pfem), m_dof(pfem)
{
    m_pMat = nullptr;
	
    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
        m_dofSU.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
        m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
    m_pMat = dynamic_cast<FEMultiphasic*>(pmat);
    assert(m_pMat);
}

//-----------------------------------------------------------------------------
// get total dof list
const FEDofList& FEMultiphasicSolidDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data.
void FEMultiphasicSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
    // get nodal DOFS
    const int nsol = m_pMat->Solutes();
    
    int N = el.Nodes();
    int ndpn = 4+nsol;
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
        
        // now the pressure dofs
        lm[ndpn*i+3] = id[m_dofP];
        
        // concentration dofs
        for (int k=0; k<nsol; ++k)
            lm[ndpn*i+4+k] = id[m_dofC+m_pMat->GetSolute(k)->GetSoluteDOF()];
        
        // rigid rotational dofs
        // TODO: Do we really need this?
        lm[ndpn*N + 3*i  ] = id[m_dofR[0]];
        lm[ndpn*N + 3*i+1] = id[m_dofR[1]];
        lm[ndpn*N + 3*i+2] = id[m_dofR[2]];
    }
    
    // substitute interface dofs for solid-shell interfaces
	FESolidElement& sel = static_cast<FESolidElement&>(el);
    for (int i=0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = m_pMesh->Node(sel.m_node[i]);
            vector<int>& id = node.m_ID;
            
            // first the back-face displacement dofs
            lm[ndpn*i  ] = id[m_dofSU[0]];
            lm[ndpn*i+1] = id[m_dofSU[1]];
            lm[ndpn*i+2] = id[m_dofSU[2]];
            
            // now the pressure dof (if the shell has it)
            if (id[m_dofQ] != -1) lm[ndpn*i+3] = id[m_dofQ];
            
            // concentration dofs
            for (int k=0; k<nsol; ++k) {
                int dofd = m_dofD+m_pMat->GetSolute(k)->GetSoluteDOF();
                if (id[dofd] != -1) lm[ndpn*i+4+k] = id[dofd];
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool FEMultiphasicSolidDomain::Init()
{
    // initialize base class
	if (FESolidDomain::Init() == false) return false;
    
    // extract the initial concentrations of the solid-bound molecules
    const int nsbm = m_pMat->SBMs();
    const int nsol = m_pMat->Solutes();
    vector<double> sbmr(nsbm, 0);
    
    for (int i = 0; i<(int)m_Elem.size(); ++i)
    {
        // get the solid element
        FESolidElement& el = m_Elem[i];
        
        // get the number of integration points
        int nint = el.GaussPoints();
        
        // loop over the integration points
        for (int n = 0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            
            // initialize sbm apparent densities
            for (int i = 0; i<nsbm; ++i)
                sbmr[i] = m_pMat->GetSBM(i)->m_rho0(mp);
            
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
            ps.m_bsb.assign(nsol, false);
            ps.m_nsbm = nsbm;
            ps.m_sbmr = sbmr;
            ps.m_sbmrp = sbmr;
            ps.m_sbmrhat.assign(nsbm,0);
            ps.m_sbmrhatp.assign(nsbm,0);
            
            // initialize referential solid volume fraction
            pb.m_phi0 = pb.m_phi0t = m_pMat->SolidReferentialVolumeFraction(mp);
            if (pb.m_phi0 > 1.0) {
                feLogError("Referential solid volume fraction of multiphasic material cannot exceed unity!\nCheck ratios of sbm apparent and true densities.");
                return false;
            }
            
            // evaluate reaction rates at initial time
            // check if this mixture includes chemical reactions
            int nreact = (int)m_pMat->Reactions();
            if (nreact) {
                // for chemical reactions involving solid-bound molecules,
                // update their concentration
                // multiphasic material point data
                FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
                
                double phi0 = pb.m_phi0t;
                for (int isbm=0; isbm<nsbm; ++isbm) {
                    // combine the molar supplies from all the reactions
                    for (int k=0; k<nreact; ++k) {
                        double zetahat = m_pMat->GetReaction(k)->ReactionSupply(mp);
                        double v = m_pMat->GetReaction(k)->m_v[nsol+isbm];
                        // remember to convert from molar supply to referential mass supply
                        ps.m_sbmrhat[isbm] += (pt.m_J-phi0)*m_pMat->SBMMolarMass(isbm)*v*zetahat;
                    }
                }
            }
        }
    }

	// set the active degrees of freedom list
	FEDofList dofs(GetFEModel());
	for (int i=0; i<nsol; ++i)
	{
		int m = m_pMat->GetSolute(i)->GetSoluteDOF();
		dofs.AddDof(m_dofC + m);
        dofs.AddDof(m_dofD + m);
	}
	m_dof = dofs;

    return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::Activate()
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
            }
        }
    }
    
    const int nsol = m_pMat->Solutes();

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
                for (int l=0; l<nsol; ++l)
                    node.set_active(m_dofD + m_pMat->GetSolute(l)->GetSoluteDOF());
            }
            else {
                node.set_active(m_dofP);
                for (int l=0; l<nsol; ++l)
                    node.set_active(m_dofC + m_pMat->GetSolute(l)->GetSoluteDOF());
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::InitMaterialPoints()
{
    const int nsol = m_pMat->Solutes();
    FEMesh& m = *GetMesh();

    // fix initial conditions for solid element nodes that are attached to the back of shells
    // this is needed because initial conditions for solid elements are prescribed to m_dofC
    // but we have to use m_dofD degrees of freedom for those solid element nodes
    for (int i=0; i<Elements(); ++i) {
        FESolidElement& el = m_Elem[i];
        // only process solid elements attached to the back of a shell
        if (el.m_bitfc.size()>0) {
            int neln = el.Nodes();
            vector<double> cic(nsol,0);
            // get the solute concentrations from nodes not attached to shells
            for (int j=0; j<neln; ++j)
            {
                FENode& node = m.Node(el.m_node[j]);
                if (!el.m_bitfc[j]) {
                    for (int l=0; l<nsol; ++l)
                        cic[l] = node.get(m_dofC + m_pMat->GetSolute(l)->GetSoluteDOF());
                    break;
                }
            }
            // assign those concentrations to nodes attached to shells
            for (int j=0; j<neln; ++j)
            {
                FENode& node = m.Node(el.m_node[j]);
                if (el.m_bitfc[j]) {
                    for (int l=0; l<nsol; ++l)
                        node.set(m_dofD + m_pMat->GetSolute(l)->GetSoluteDOF(), cic[l]);
                }
            }
        }
    }

    const int nsbm = m_pMat->SBMs();
    
    const int NE = FEElement::MAX_NODES;
    double p0[NE];
    vector< vector<double> > c0(nsol, vector<double>(NE));
    vector<int> sid(nsol);
    for (int j = 0; j<nsol; ++j) sid[j] = m_pMat->GetSolute(j)->GetSoluteDOF();
    
    for (int j = 0; j<(int)m_Elem.size(); ++j)
    {
        // get the solid element
        FESolidElement& el = m_Elem[j];
        
        // get the number of nodes
        int neln = el.Nodes();
        // get initial values of fluid pressure and solute concentrations
        if (el.m_bitfc.size() == 0) {
            for (int i = 0; i<neln; ++i)
            {
                FENode& ni = m.Node(el.m_node[i]);
                p0[i] = ni.get(m_dofP);
                for (int isol = 0; isol<nsol; ++isol)
                    c0[isol][i] = ni.get(m_dofC + sid[isol]);
            }
        }
        else {
            for (int i = 0; i<neln; ++i)
            {
                FENode& ni = m.Node(el.m_node[i]);
                p0[i] = el.m_bitfc[i] ? ni.get(m_dofQ) : ni.get(m_dofP);
                for (int isol = 0; isol<nsol; ++isol)
                    c0[isol][i] = (ni.m_ID[m_dofD + isol] != -1) ? ni.get(m_dofD + sid[isol]) : ni.get(m_dofC + sid[isol]);
            }
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
            
            // initialize multiphasic solutes
            ps.m_nsol = nsol;
            ps.m_nsbm = nsbm;
            
            // initialize effective solute concentrations
            for (int isol = 0; isol<nsol; ++isol) {
                ps.m_c[isol] = el.Evaluate(c0[isol], n);
                ps.m_gradc[isol] = gradient(el, c0[isol], n);
            }
            
            ps.m_psi = m_pMat->ElectricPotential(mp);
            for (int isol = 0; isol<nsol; ++isol) {
                ps.m_ca[isol] = m_pMat->Concentration(mp, isol);
                ps.m_j[isol] = m_pMat->SoluteFlux(mp, isol);
                ps.m_crp[isol] = pm.m_J*m_pMat->Porosity(mp)*ps.m_ca[isol];
            }
            pt.m_pa = m_pMat->Pressure(mp);
            
            // determine if solute is 'solid-bound'
            for (int isol = 0; isol<nsol; ++isol) {
                FESolute* soli = m_pMat->GetSolute(isol);
                if (soli->m_pDiff->Diffusivity(mp).norm() == 0) ps.m_bsb[isol] = true;
            }
            
            // initialize referential solid volume fraction
            pt.m_phi0t = m_pMat->SolidReferentialVolumeFraction(mp);
            
            // calculate FCD, current and stress
            ps.m_cF = m_pMat->FixedChargeDensity(mp);
            ps.m_Ie = m_pMat->CurrentDensity(mp);
            pm.m_s = m_pMat->Stress(mp);
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::Reset()
{
    // reset base class
    FESolidDomain::Reset();
    
    const int nsol = m_pMat->Solutes();
    const int nsbm = m_pMat->SBMs();
    
    // extract the initial concentrations of the solid-bound molecules
    vector<double> sbmr(nsbm,0);
    
    for (int i=0; i<(int) m_Elem.size(); ++i)
    {
        // get the solid element
        FESolidElement& el = m_Elem[i];
        
        // get the number of integration points
        int nint = el.GaussPoints();
        
        // loop over the integration points
        for (int n=0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            
            // initialize sbm apparent densities
            for (int i = 0; i<nsbm; ++i)
                sbmr[i] = m_pMat->GetSBM(i)->m_rho0(mp);
            
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
            ps.m_bsb.assign(nsol, false);
            ps.m_nsbm = nsbm;
            ps.m_sbmr = sbmr;
            ps.m_sbmrp = sbmr;
            ps.m_sbmrhat.assign(nsbm,0);
            ps.m_sbmrhatp.assign(nsbm,0);

            // initialize referential solid volume fraction
            pt.m_phi0 = pt.m_phi0t = m_pMat->SolidReferentialVolumeFraction(mp);

            // reset chemical reaction element data
            ps.m_cri.clear();
            ps.m_crd.clear();
            for (int j=0; j<m_pMat->Reactions(); ++j)
                m_pMat->GetReaction(j)->ResetElementData(mp);
        }
    }
    m_breset = true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    FESolidDomain::PreSolveUpdate(timeInfo);
    
    const int NE = FEElement::MAX_NODES;
    vec3d x0[NE], xt[NE], r0, rt;
    FEMesh& m = *GetMesh();
    for (size_t iel=0; iel<m_Elem.size(); ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        int neln = el.Nodes();
        for (int i=0; i<neln; ++i)
        {
            x0[i] = m.Node(el.m_node[i]).m_r0;
            xt[i] = m.Node(el.m_node[i]).m_rt;
        }
        
        int n = el.GaussPoints();
        for (int j=0; j<n; ++j)
        {
            r0 = el.Evaluate(x0, j);
            rt = el.Evaluate(xt, j);
            
            FEMaterialPoint& mp = *el.GetMaterialPoint(j);
            FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
            FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            FEMultigenSBMMaterialPoint* pmg = mp.ExtractData<FEMultigenSBMMaterialPoint>();
            
            mp.m_r0 = r0;
            mp.m_rt = rt;
            
            pe.m_J = defgrad(el, pe.m_F, j);
            
            // reset determinant of solid deformation gradient at previous time
            pt.m_Jp = pe.m_J;
            
            // reset referential solid volume fraction at previous time
            pt.m_phi0p = pt.m_phi0t;
            
            // reset referential actual solute concentration at previous time
            for (int j=0; j<m_pMat->Solutes(); ++j) {
                ps.m_crp[j] = pe.m_J*m_pMat->Porosity(mp)*ps.m_ca[j];
            }
            
            // reset referential solid-bound molecule concentrations at previous time
            ps.m_sbmrp = ps.m_sbmr;
            ps.m_sbmrhatp = ps.m_sbmrhat;
            
            // reset generational referential solid-bound molecule concentrations at previous time
            if (pmg) {
                for (int i=0; i<pmg->m_ngen; ++i) {
                    for (int j=0; j<ps.m_nsbm; ++j) {
                        pmg->m_gsbmrp[i][j] = pmg->m_gsbmr[i][j];
                    }
                }
            }
            
            // reset chemical reaction element data
            for (int j=0; j<m_pMat->Reactions(); ++j)
                m_pMat->GetReaction(j)->InitializeElementData(mp);
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::InternalForces(FEGlobalVector& R)
{
    size_t NE = m_Elem.size();
    
    // get nodal DOFS
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = ndpn*el.Nodes();
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

void FEMultiphasicSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
    int i, isol, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradN;
    mat3ds s;
    
    const double* Gr, *Gs, *Gt, *H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*	gw = el.GaussWeights();
    
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    const int nreact = m_pMat->Reactions();
    
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
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // next we get the determinant
        double Jp = bpt.m_Jp;
        double J = pt.m_J;
        
        // and then finally
        double divv = ((J-Jp)/dt)/J;
        
        // get the stress for this integration point
        s = pt.m_s;
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        vector<vec3d> j(spt.m_j);
        vector<int> z(nsol);
        vector<double> kappa(spt.m_k);
        vec3d je(0,0,0);
        
        for (isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            je += j[isol]*z[isol];
        }
        
        // evaluate the porosity, its derivative w.r.t. J, and its gradient
        double phiw = m_pMat->Porosity(mp);
        vector<double> chat(nsol,0);
        
        // get the solvent supply
        double phiwhat = 0;
        if (m_pMat->GetSolventSupply()) phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
        
        // chemical reactions
        for (i=0; i<nreact; ++i) {
            FEChemicalReaction* pri = m_pMat->GetReaction(i);
            double zhat = pri->ReactionSupply(mp);
            phiwhat += phiw*pri->m_Vbar*zhat;
            for (isol=0; isol<nsol; ++isol)
                chat[isol] += phiw*zhat*pri->m_v[isol];
        }
        
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
            fe[ndpn*i  ] -= fu.x*detJt;
            fe[ndpn*i+1] -= fu.y*detJt;
            fe[ndpn*i+2] -= fu.z*detJt;
            fe[ndpn*i+3] -= dt*(w*gradN + (phiwhat - divv)*H[i])*detJt;
            for (isol=0; isol<nsol; ++isol)
                fe[ndpn*i+4+isol] -= dt*(gradN*(j[isol]+je*m_pMat->m_penalty)
                                         + H[i]*(chat[isol] - (phiw*spt.m_ca[isol] - spt.m_crp[isol]/J)/dt)
                                         )*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::InternalForcesSS(FEGlobalVector& R)
{
    size_t NE = m_Elem.size();
    
    // get nodal DOFS
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = ndpn*el.Nodes();
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

void FEMultiphasicSolidDomain::ElementInternalForceSS(FESolidElement& el, vector<double>& fe)
{
    int i, isol, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradN;
    mat3ds s;
    
    const double* Gr, *Gs, *Gt, *H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*	gw = el.GaussWeights();
    
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    const int nreact = m_pMat->Reactions();
    
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
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // get the stress for this integration point
        s = pt.m_s;
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        vector<vec3d> j(spt.m_j);
        vector<int> z(nsol);
        vector<double> kappa(spt.m_k);
        vec3d je(0,0,0);
        
        for (isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            je += j[isol]*z[isol];
        }
        
        // evaluate the porosity, its derivative w.r.t. J, and its gradient
        double phiw = m_pMat->Porosity(mp);
        vector<double> chat(nsol,0);
        
        // get the solvent supply
        double phiwhat = 0;
        if (m_pMat->GetSolventSupply()) phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
        
        // chemical reactions
        for (i=0; i<nreact; ++i) {
            FEChemicalReaction* pri = m_pMat->GetReaction(i);
            double zhat = pri->ReactionSupply(mp);
            phiwhat += phiw*pri->m_Vbar*zhat;
            for (isol=0; isol<nsol; ++isol)
                chat[isol] += phiw*zhat*pri->m_v[isol];
        }
        
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
            fe[ndpn*i  ] -= fu.x*detJt;
            fe[ndpn*i+1] -= fu.y*detJt;
            fe[ndpn*i+2] -= fu.z*detJt;
            fe[ndpn*i+3] -= dt*(w*gradN + H[i]*phiwhat)*detJt;
            for (isol=0; isol<nsol; ++isol)
                fe[ndpn*i+4+isol] -= dt*(gradN*(j[isol]+je*m_pMat->m_penalty)
                                         + H[i]*phiw*chat[isol]
                                         )*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::StiffnessMatrix(FELinearSystem& LS, bool bsymm)
{
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FESolidElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);

        // allocate stiffness matrix
        int neln = el.Nodes();
        int ndof = neln*ndpn;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementMultiphasicStiffness(el, ke, bsymm);

		// get the lm vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::StiffnessMatrixSS(FELinearSystem& LS, bool bsymm)
{
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FESolidElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);

        // allocate stiffness matrix
        int neln = el.Nodes();
        int ndof = neln*ndpn;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementMultiphasicStiffnessSS(el, ke, bsymm);

		// get the lm vector
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
bool FEMultiphasicSolidDomain::ElementMultiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm)
{
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(neln);
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    const int nsbm   = m_pMat->SBMs();
    const int nreact = m_pMat->Reactions();
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (int n=0; n<nint; ++n)
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
        for (int i=0; i<neln; ++i)
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
        
        vector<double> c(spt.m_c);
        vector<vec3d> gradc(spt.m_gradc);
        vector<int> z(nsol);
        
        vector<double> kappa(spt.m_k);
        
        // get the charge number
        for (int isol=0; isol<nsol; ++isol)
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
        
        vector<double> dkdJ(spt.m_dkdJ);
        vector< vector<double> > dkdc(spt.m_dkdc);
        vector< vector<double> > dkdr(spt.m_dkdr);
        vector< vector<double> > dkdJr(spt.m_dkdJr);
        vector< vector< vector<double> > > dkdrc(spt.m_dkdrc);
        
        // evaluate the porosity and its derivative
        double phiw = m_pMat->Porosity(mp);
        double phi0 = ppt.m_phi0t;
        double phis = 1. - phiw;
        double dpdJ = phis/J;
        
        // evaluate the osmotic coefficient
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        
        // evaluate the permeability
        mat3ds K = m_pMat->GetPermeability()->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        vector<mat3ds> dKdc(nsol);
        vector<mat3ds> D(nsol);
        vector<tens4dmm> dDdE(nsol);
        vector< vector<mat3ds> > dDdc(nsol, vector<mat3ds>(nsol));
        vector<double> D0(nsol);
        vector< vector<double> > dD0dc(nsol, vector<double>(nsol));
        vector<double> dodc(nsol);
        vector<mat3ds> dTdc(nsol);
        vector<mat3ds> ImD(nsol);
        mat3dd I(1);
        
        // evaluate the solvent supply and its derivatives
        mat3ds Phie; Phie.zero();
        double Phip = 0;
        vector<double> Phic(nsol,0);
        vector<mat3ds> dchatde(nsol);
        if (m_pMat->GetSolventSupply()) {
            Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
            Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
        }
        
        // chemical reactions
		vector<double> reactionSupply(nreact, 0.0);
		vector<mat3ds> tangentReactionSupplyStrain(nreact);
		vector< vector<double> > tangentReactionSupplyConcentration(nreact, vector<double>(nsol));
		for (int i = 0; i < nreact; ++i)
		{
			FEChemicalReaction* reacti = m_pMat->GetReaction(i);

			reactionSupply[i] = reacti->ReactionSupply(mp);
			tangentReactionSupplyStrain[i] = reacti->Tangent_ReactionSupply_Strain(mp);

			for (int isol = 0; isol < nsol; ++isol)
			{
				tangentReactionSupplyConcentration[i][isol] = reacti->Tangent_ReactionSupply_Concentration(mp, isol);
			}

			Phie += reacti->m_Vbar*(I*reactionSupply[i]
				+ tangentReactionSupplyStrain [i]*(J*phiw));

		}
        
        for (int isol=0; isol<nsol; ++isol) {
        
			FESolute* soli = m_pMat->GetSolute(isol);

			// evaluate the permeability derivatives
            dKdc[isol] = m_pMat->GetPermeability()->Tangent_Permeability_Concentration(mp,isol);
            
            // evaluate the diffusivity tensor and its derivatives
            D[isol] = soli->m_pDiff->Diffusivity(mp);
            dDdE[isol] = soli->m_pDiff->Tangent_Diffusivity_Strain(mp);
            
            // evaluate the solute free diffusivity
            D0[isol] = soli->m_pDiff->Free_Diffusivity(mp);
            
            // evaluate the derivative of the osmotic coefficient
            dodc[isol] = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp,isol);
            
            // evaluate the stress tangent with concentration
            //			dTdc[isol] = pm->GetSolid()->Tangent_Concentration(mp,isol);
            dTdc[isol] = mat3ds(0,0,0,0,0,0);
            
            ImD[isol] = I-D[isol]/D0[isol];
            
            for (int jsol=0; jsol<nsol; ++jsol) {
                dDdc[isol][jsol] = soli->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
                dD0dc[isol][jsol] = soli->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
            }
            
            // evaluate the solvent supply tangent with concentration
            if (m_pMat->GetSolventSupply()) Phic[isol] = m_pMat->GetSolventSupply()->Tangent_Supply_Concentration(mp,isol);
            
            // chemical reactions
            dchatde[isol].zero();
            for (int ireact=0; ireact<nreact; ++ireact) {
				FEChemicalReaction* reacti = m_pMat->GetReaction(ireact);

                dchatde[isol] += reacti->m_v[isol]
					*(I*reactionSupply[ireact]
					+ tangentReactionSupplyStrain[ireact] *(J*phiw));

                Phic[isol] += phiw* reacti->m_Vbar*tangentReactionSupplyConcentration[ireact][isol];
            }
        }
        
        // Miscellaneous constants
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double penalty = m_pMat->m_penalty;
        
        // evaluate the effective permeability and its derivatives
        mat3ds Ki = K.inverse();
        mat3ds Ke(0,0,0,0,0,0);
        tens4d G = (dyad1(Ki,I) - dyad4(Ki,I)*2)*2 - ddot(dyad2(Ki,Ki),dKdE);
        vector<mat3ds> Gc(nsol);
        vector<mat3ds> dKedc(nsol);
        for (int isol=0; isol<nsol; ++isol) {
            Ke += ImD[isol]*(kappa[isol]*c[isol]/D0[isol]);
            G += dyad1(ImD[isol],I)*(R*T*c[isol]*J/D0[isol]/phiw*(dkdJ[isol]-kappa[isol]/phiw*dpdJ))
            +(dyad1(I,I) - dyad2(I,I)*2 - dDdE[isol]/D0[isol])*(R*T*kappa[isol]*c[isol]/phiw/D0[isol]);
            Gc[isol] = ImD[isol]*(kappa[isol]/D0[isol]);
            for (int jsol=0; jsol<nsol; ++jsol) {
                Gc[isol] += ImD[jsol]*(c[jsol]/D0[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol]))
                -(dDdc[jsol][isol]-D[jsol]*(dD0dc[jsol][isol]/D0[jsol])*(kappa[jsol]*c[jsol]/SQR(D0[jsol])));
            }
            Gc[isol] *= R*T/phiw;
        }
        Ke = (Ki + Ke*(R*T/phiw)).inverse();
        tens4d dKedE = (dyad1(Ke,I) - dyad4(Ke,I)*2)*2 - ddot(dyad2(Ke,Ke),G);
        for (int isol=0; isol<nsol; ++isol)
            dKedc[isol] = -(Ke*(-Ki*dKdc[isol]*Ki + Gc[isol])*Ke).sym();
        
        // calculate all the matrices
        vec3d vtmp,gp,qpu;
        vector<vec3d> gc(nsol),qcu(nsol),wc(nsol),jce(nsol);
        vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
        mat3d wu, jue;
        vector<mat3d> ju(nsol);
        vector< vector<double> > qcc(nsol, vector<double>(nsol));
        vector< vector<double> > dchatdc(nsol, vector<double>(nsol));
        double sum;
        mat3ds De;
        for (int i=0; i<neln; ++i)
        {
            for (int j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradN[i]*(s*gradN[j])) + vdotTdotv(gradN[i], C, gradN[j]))*detJ;
                ke[ndpn*i  ][ndpn*j  ] += Kuu[0][0]; ke[ndpn*i  ][ndpn*j+1] += Kuu[0][1]; ke[ndpn*i  ][ndpn*j+2] += Kuu[0][2];
                ke[ndpn*i+1][ndpn*j  ] += Kuu[1][0]; ke[ndpn*i+1][ndpn*j+1] += Kuu[1][1]; ke[ndpn*i+1][ndpn*j+2] += Kuu[1][2];
                ke[ndpn*i+2][ndpn*j  ] += Kuu[2][0]; ke[ndpn*i+2][ndpn*j+1] += Kuu[2][1]; ke[ndpn*i+2][ndpn*j+2] += Kuu[2][2];
                
                // calculate the kpu matrix
                gp = vec3d(0,0,0);
                for (int isol=0; isol<nsol; ++isol) gp += (D[isol]*gradc[isol])*(kappa[isol]/D0[isol]);
                gp = gradp+gp*(R*T);
                wu = vdotTdotv(-gp, dKedE, gradN[j]);
                for (int isol=0; isol<nsol; ++isol) {
                    wu += (((Ke*(D[isol]*gradc[isol])) & gradN[j])*(J*dkdJ[isol] - kappa[isol])
                           +Ke*(2*kappa[isol]*(gradN[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
                    + (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradN[j]))*(-kappa[isol]*R*T/D0[isol]);
                }
                qpu = -gradN[j]*(1.0/dt);
                vtmp = (wu.transpose()*gradN[i] + (qpu + Phie*gradN[j])*H[i])*(detJ*dt);
                ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
                ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
                ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
                
                // calculate the kup matrix
                vtmp = -gradN[i]*H[j]*detJ;
                ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
                ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
                ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
                
                // calculate the kpp matrix
                ke[ndpn*i+3][ndpn*j+3] += (H[i]*H[j]*Phip - gradN[i]*(Ke*gradN[j]))*(detJ*dt);
                
                // calculate kcu matrix data
                jue.zero();
                De.zero();
                for (int isol=0; isol<nsol; ++isol) {
                    gc[isol] = -gradc[isol]*phiw + w*c[isol]/D0[isol];
                    ju[isol] = ((D[isol]*gc[isol]) & gradN[j])*(J*dkdJ[isol])
                    + vdotTdotv(gc[isol], dDdE[isol], gradN[j])*kappa[isol]
                    + (((D[isol]*gradc[isol]) & gradN[j])*(-phis)
                       +(D[isol]*((gradN[j]*w)*2) - ((D[isol]*w) & gradN[j]))*c[isol]/D0[isol]
                       )*kappa[isol]
                    +D[isol]*wu*(kappa[isol]*c[isol]/D0[isol]);
                    jue += ju[isol]*z[isol];
                    De += D[isol]*(z[isol]*kappa[isol]*c[isol]/D0[isol]);
                    qcu[isol] = qpu*(c[isol]*(kappa[isol]+J*phiw*dkdJ[isol]));
                    
                    // chemical reactions
                    for (int ireact=0; ireact<nreact; ++ireact) {
						FEChemicalReaction* reacti = m_pMat->GetReaction(ireact);

                        double sum1 = 0;
                        double sum2 = 0;
                        for (int isbm=0; isbm<nsbm; ++isbm) {
                            sum1 += m_pMat->SBMMolarMass(isbm)*reacti->m_v[nsol+isbm]*
                            ((J-phi0)*dkdr[isol][isbm]-kappa[isol]/m_pMat->SBMDensity(isbm));
                            sum2 += m_pMat->SBMMolarMass(isbm)*reacti->m_v[nsol+isbm]*
                            (dkdr[isol][isbm]+(J-phi0)*dkdJr[isol][isbm]-dkdJ[isol]/m_pMat->SBMDensity(isbm));
                        }
                        double zhat = reactionSupply[ireact];
                        mat3dd zhatI(zhat);
                        mat3ds dzde = tangentReactionSupplyStrain[ireact];
                        qcu[isol] -= ((zhatI+dzde*(J-phi0))*gradN[j])*(sum1*c[isol])
                        +gradN[j]*(c[isol]*(J-phi0)*sum2*zhat);
                    }
                }
                
                for (int isol=0; isol<nsol; ++isol) {
                    
                    // calculate the kcu matrix
                    vtmp = ((ju[isol]+jue*penalty).transpose()*gradN[i]
                            + (qcu[isol] + dchatde[isol]*gradN[j])*H[i])*(detJ*dt);
                    ke[ndpn*i+4+isol][ndpn*j  ] += vtmp.x;
                    ke[ndpn*i+4+isol][ndpn*j+1] += vtmp.y;
                    ke[ndpn*i+4+isol][ndpn*j+2] += vtmp.z;
                    
                    // calculate the kcp matrix
                    ke[ndpn*i+4+isol][ndpn*j+3] -= (gradN[i]*(
                                                              (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                               +De*penalty)
                                                              *(Ke*gradN[j])
                                                              ))*(detJ*dt);
                    
                    // calculate the kuc matrix
                    sum = 0;
                    for (int jsol=0; jsol<nsol; ++jsol)
                        sum += c[jsol]*(dodc[isol]*kappa[jsol]+osmc*dkdc[jsol][isol]);
                    vtmp = (dTdc[isol]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[isol]+sum)))*H[j]*detJ;
                    ke[ndpn*i  ][ndpn*j+4+isol] += vtmp.x;
                    ke[ndpn*i+1][ndpn*j+4+isol] += vtmp.y;
                    ke[ndpn*i+2][ndpn*j+4+isol] += vtmp.z;
                    
                    // calculate the kpc matrix
                    vtmp = vec3d(0,0,0);
                    for (int jsol=0; jsol<nsol; ++jsol)
                        vtmp += (D[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol])
                                 +dDdc[jsol][isol]*kappa[jsol])/D0[jsol]*gradc[jsol];
                    wc[isol] = (dKedc[isol]*gp)*(-H[j])
                    -Ke*((D[isol]*gradN[j])*(kappa[isol]/D0[isol])+vtmp*H[j])*(R*T);
                    ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol]+H[i]*H[j]*Phic[isol])*(detJ*dt);
                    
                }
                
                // calculate data for the kcc matrix
                jce.assign(nsol, vec3d(0,0,0));
                for (int isol=0; isol<nsol; ++isol) {
                    for (int jsol=0; jsol<nsol; ++jsol) {
                        if (jsol != isol) {
                            jc[isol][jsol] =
                            ((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            
                            qcc[isol][jsol] = -H[j]*phiw/dt*c[isol]*dkdc[isol][jsol];
                        }
                        else {
                            jc[isol][jsol] = (D[isol]*(gradN[j]*(-phiw)+w*(H[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            
                            qcc[isol][jsol] = -H[j]*phiw/dt*(c[isol]*dkdc[isol][jsol] + kappa[isol]);
                        }
                        jce[jsol] += jc[isol][jsol]*z[isol];
                        
                        // chemical reactions
                        dchatdc[isol][jsol] = 0;
                        for (int ireact=0; ireact<nreact; ++ireact) {
							FEChemicalReaction* reacti = m_pMat->GetReaction(ireact);

                            dchatdc[isol][jsol] += reacti->m_v[isol]
                            * tangentReactionSupplyConcentration[ireact][jsol];

                            double sum1 = 0;
                            double sum2 = 0;
                            for (int isbm=0; isbm<nsbm; ++isbm) {
                                sum1 += m_pMat->SBMMolarMass(isbm)*reacti->m_v[nsol+isbm]*
                                ((J-phi0)*dkdr[isol][isbm]-kappa[isol]/m_pMat->SBMDensity(isbm));
                                sum2 += m_pMat->SBMMolarMass(isbm)*reacti->m_v[nsol+isbm]*
                                ((J-phi0)*dkdrc[isol][isbm][jsol]-dkdc[isol][jsol]/m_pMat->SBMDensity(isbm));
                            }
                            double zhat = reactionSupply[ireact];
                            double dzdc = tangentReactionSupplyConcentration[ireact][jsol];
                            if (jsol != isol) {
                                qcc[isol][jsol] -= H[j]*phiw*c[isol]*(dzdc*sum1+zhat*sum2);
                            }
                            else {
                                qcc[isol][jsol] -= H[j]*phiw*((zhat+c[isol]*dzdc)*sum1+c[isol]*zhat*sum2);
                            }
                        }
                    }
                }
                
                // calculate the kcc matrix
                for (int isol=0; isol<nsol; ++isol) {
                    for (int jsol=0; jsol<nsol; ++jsol) {
                        ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty)
                                                             + H[i]*(qcc[isol][jsol]
                                                                     + H[j]*phiw*dchatdc[isol][jsol]))*(detJ*dt);
                    }
                }
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
    double tmp;
    if (bsymm) {
        for (int i=0; i<ndpn*neln; ++i)
            for (int j=i+1; j<ndpn*neln; ++j) {
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
bool FEMultiphasicSolidDomain::ElementMultiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm)
{
    int i, j, isol, jsol, n, ireact;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(neln);
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    const int nreact = m_pMat->Reactions();
    
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
        
        vector<double> c(spt.m_c);
        vector<vec3d> gradc(spt.m_gradc);
        vector<int> z(nsol);
        
        vector<double> zz(nsol);
        vector<double> kappa(spt.m_k);
        
        // get the charge number
        for (isol=0; isol<nsol; ++isol)
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
        
        vector<double> dkdJ(spt.m_dkdJ);
        vector< vector<double> > dkdc(spt.m_dkdc);
        
        // evaluate the porosity and its derivative
        double phiw = m_pMat->Porosity(mp);
        double phis = 1. - phiw;
        double dpdJ = phis/J;
        
        // evaluate the osmotic coefficient
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        
        // evaluate the permeability
        mat3ds K = m_pMat->GetPermeability()->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        vector<mat3ds> dKdc(nsol);
        vector<mat3ds> D(nsol);
        vector<tens4dmm> dDdE(nsol);
        vector< vector<mat3ds> > dDdc(nsol, vector<mat3ds>(nsol));
        vector<double> D0(nsol);
        vector< vector<double> > dD0dc(nsol, vector<double>(nsol));
        vector<double> dodc(nsol);
        vector<mat3ds> dTdc(nsol);
        vector<mat3ds> ImD(nsol);
        mat3dd I(1);
        
        // evaluate the solvent supply and its derivatives
        double phiwhat = 0;
        mat3ds Phie; Phie.zero();
        double Phip = 0;
        vector<double> Phic(nsol,0);
        if (m_pMat->GetSolventSupply()) {
            phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
            Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
            Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
        }
        
        // chemical reactions
        for (i=0; i<nreact; ++i)
            Phie += m_pMat->GetReaction(i)->m_Vbar*(I*m_pMat->GetReaction(i)->ReactionSupply(mp)
                                                    +m_pMat->GetReaction(i)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
        
        for (isol=0; isol<nsol; ++isol) {
            // evaluate the permeability derivatives
            dKdc[isol] = m_pMat->GetPermeability()->Tangent_Permeability_Concentration(mp,isol);
            
            // evaluate the diffusivity tensor and its derivatives
            D[isol] = m_pMat->GetSolute(isol)->m_pDiff->Diffusivity(mp);
            dDdE[isol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Diffusivity_Strain(mp);
            
            // evaluate the solute free diffusivity
            D0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            
            // evaluate the derivative of the osmotic coefficient
            dodc[isol] = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp,isol);
            
            // evaluate the stress tangent with concentration
            //			dTdc[isol] = pm->GetSolid()->Tangent_Concentration(mp,isol);
            dTdc[isol] = mat3ds(0,0,0,0,0,0);
            
            ImD[isol] = I-D[isol]/D0[isol];
            
            for (jsol=0; jsol<nsol; ++jsol) {
                dDdc[isol][jsol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
                dD0dc[isol][jsol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
            }
            
            // evaluate the solvent supply tangent with concentration
            if (m_pMat->GetSolventSupply()) Phic[isol] = m_pMat->GetSolventSupply()->Tangent_Supply_Concentration(mp,isol);
            
        }
        
        // Miscellaneous constants
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double penalty = m_pMat->m_penalty;
        
        // evaluate the effective permeability and its derivatives
        mat3ds Ki = K.inverse();
        mat3ds Ke(0,0,0,0,0,0);
        tens4d G = (dyad1(Ki,I) - dyad4(Ki,I)*2)*2 - ddot(dyad2(Ki,Ki),dKdE);
        vector<mat3ds> Gc(nsol);
        vector<mat3ds> dKedc(nsol);
        for (isol=0; isol<nsol; ++isol) {
            Ke += ImD[isol]*(kappa[isol]*c[isol]/D0[isol]);
            G += dyad1(ImD[isol],I)*(R*T*c[isol]*J/D0[isol]/phiw*(dkdJ[isol]-kappa[isol]/phiw*dpdJ))
            +(dyad1(I,I) - dyad2(I,I)*2 - dDdE[isol]/D0[isol])*(R*T*kappa[isol]*c[isol]/phiw/D0[isol]);
            Gc[isol] = ImD[isol]*(kappa[isol]/D0[isol]);
            for (jsol=0; jsol<nsol; ++jsol) {
                Gc[isol] += ImD[jsol]*(c[jsol]/D0[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol]))
                -(dDdc[jsol][isol]-D[jsol]*(dD0dc[jsol][isol]/D0[jsol])*(kappa[jsol]*c[jsol]/SQR(D0[jsol])));
            }
            Gc[isol] *= R*T/phiw;
        }
        Ke = (Ki + Ke*(R*T/phiw)).inverse();
        tens4d dKedE = (dyad1(Ke,I) - 2*dyad4(Ke,I))*2 - ddot(dyad2(Ke,Ke),G);
        for (isol=0; isol<nsol; ++isol)
            dKedc[isol] = -(Ke*(-Ki*dKdc[isol]*Ki + Gc[isol])*Ke).sym();
        
        // calculate all the matrices
        vec3d vtmp,gp,qpu;
        vector<vec3d> gc(nsol),wc(nsol),jce(nsol);
        vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
        mat3d wu, jue;
        vector<mat3d> ju(nsol);
        vector< vector<double> > dchatdc(nsol, vector<double>(nsol));
        double sum;
        mat3ds De;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradN[i]*(s*gradN[j])) + vdotTdotv(gradN[i], C, gradN[j]))*detJ;
                ke[ndpn*i  ][ndpn*j  ] += Kuu[0][0]; ke[ndpn*i  ][ndpn*j+1] += Kuu[0][1]; ke[ndpn*i  ][ndpn*j+2] += Kuu[0][2];
                ke[ndpn*i+1][ndpn*j  ] += Kuu[1][0]; ke[ndpn*i+1][ndpn*j+1] += Kuu[1][1]; ke[ndpn*i+1][ndpn*j+2] += Kuu[1][2];
                ke[ndpn*i+2][ndpn*j  ] += Kuu[2][0]; ke[ndpn*i+2][ndpn*j+1] += Kuu[2][1]; ke[ndpn*i+2][ndpn*j+2] += Kuu[2][2];
                
                // calculate the kpu matrix
                gp = vec3d(0,0,0);
                for (isol=0; isol<nsol; ++isol) gp += (D[isol]*gradc[isol])*(kappa[isol]/D0[isol]);
                gp = gradp+gp*(R*T);
                wu = vdotTdotv(-gp, dKedE, gradN[j]);
                for (isol=0; isol<nsol; ++isol) {
                    wu += (((Ke*(D[isol]*gradc[isol])) & gradN[j])*(J*dkdJ[isol] - kappa[isol])
                           +Ke*(2*kappa[isol]*(gradN[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
                    + (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradN[j]))*(-kappa[isol]*R*T/D0[isol]);
                }
                qpu = Phie*gradN[j];
                vtmp = (wu.transpose()*gradN[i] + qpu*H[i])*(detJ*dt);
                ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
                ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
                ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
                
                // calculate the kup matrix
                vtmp = -gradN[i]*H[j]*detJ;
                ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
                ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
                ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
                
                // calculate the kpp matrix
                ke[ndpn*i+3][ndpn*j+3] += (H[i]*H[j]*Phip - gradN[i]*(Ke*gradN[j]))*(detJ*dt);
                
                // calculate kcu matrix data
                jue.zero();
                De.zero();
                for (isol=0; isol<nsol; ++isol) {
                    gc[isol] = -gradc[isol]*phiw + w*c[isol]/D0[isol];
                    ju[isol] = ((D[isol]*gc[isol]) & gradN[j])*(J*dkdJ[isol])
                    + vdotTdotv(gc[isol], dDdE[isol], gradN[j])*kappa[isol]
                    + (((D[isol]*gradc[isol]) & gradN[j])*(-phis)
                       +(D[isol]*((gradN[j]*w)*2) - ((D[isol]*w) & gradN[j]))*c[isol]/D0[isol]
                       )*kappa[isol]
                    +D[isol]*wu*(kappa[isol]*c[isol]/D0[isol]);
                    jue += ju[isol]*z[isol];
                    De += D[isol]*(z[isol]*kappa[isol]*c[isol]/D0[isol]);
                }
                
                for (isol=0; isol<nsol; ++isol) {
                    
                    // calculate the kcu matrix
                    vtmp = ((ju[isol]+jue*penalty).transpose()*gradN[i])*(detJ*dt);
                    ke[ndpn*i+4+isol][ndpn*j  ] += vtmp.x;
                    ke[ndpn*i+4+isol][ndpn*j+1] += vtmp.y;
                    ke[ndpn*i+4+isol][ndpn*j+2] += vtmp.z;
                    
                    // calculate the kcp matrix
                    ke[ndpn*i+4+isol][ndpn*j+3] -= (gradN[i]*(
                                                              (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                               +De*penalty)
                                                              *(Ke*gradN[j])
                                                              ))*(detJ*dt);
                    
                    // calculate the kuc matrix
                    sum = 0;
                    for (jsol=0; jsol<nsol; ++jsol)
                        sum += c[jsol]*(dodc[isol]*kappa[jsol]+osmc*dkdc[jsol][isol]);
                    vtmp = (dTdc[isol]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[isol]+sum)))*H[j]*detJ;
                    ke[ndpn*i  ][ndpn*j+4+isol] += vtmp.x;
                    ke[ndpn*i+1][ndpn*j+4+isol] += vtmp.y;
                    ke[ndpn*i+2][ndpn*j+4+isol] += vtmp.z;
                    
                    // calculate the kpc matrix
                    vtmp = vec3d(0,0,0);
                    for (jsol=0; jsol<nsol; ++jsol)
                        vtmp += (D[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol])
                                 +dDdc[jsol][isol]*kappa[jsol])/D0[jsol]*gradc[jsol];
                    wc[isol] = (dKedc[isol]*gp)*(-H[j])
                    -Ke*((D[isol]*gradN[j])*(kappa[isol]/D0[isol])+vtmp*H[j])*(R*T);
                    ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol]+H[i]*H[j]*Phic[isol])*(detJ*dt);
                    
                }
                
                // calculate data for the kcc matrix
                jce.assign(nsol, vec3d(0,0,0));
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        if (jsol != isol) {
                            jc[isol][jsol] = 
                            ((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                        }
                        else {
                            jc[isol][jsol] = (D[isol]*(gradN[j]*(-phiw)+w*(H[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                        }
                        jce[jsol] += jc[isol][jsol]*z[isol];
                        
                        // chemical reactions
                        dchatdc[isol][jsol] = 0;
                        for (ireact=0; ireact<nreact; ++ireact)
                            dchatdc[isol][jsol] += m_pMat->GetReaction(ireact)->m_v[isol]
                            *m_pMat->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp,jsol);
                    }
                }
                
                // calculate the kcc matrix
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty)
                                                             + H[i]*H[j]*phiw*dchatdc[isol][jsol])*(detJ*dt);
                    }
                }
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
    double tmp;
    if (bsymm) {
        for (i=0; i<ndpn*neln; ++i)
            for (j=i+1; j<ndpn*neln; ++j) {
                tmp = 0.5*(ke[i][j]+ke[j][i]);
                ke[i][j] = ke[j][i] = tmp;
            }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicSolidDomain::Update(const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();
    bool berr = false;
    int NE = (int) m_Elem.size();
    double dt = fem.GetTime().timeIncrement;
#pragma omp parallel for shared(NE, berr)
    for (int i=0; i<NE; ++i)
    {
        try
        {
            UpdateElementStress(i, dt);
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
void FEMultiphasicSolidDomain::UpdateElementStress(int iel, double dt)
{
    int j, k, n;
    int nint, neln;
    double* gw;
    vec3d r0[FEElement::MAX_NODES];
    vec3d rt[FEElement::MAX_NODES];
    double pn[FEElement::MAX_NODES];
    
    FEMesh& mesh = *m_pMesh;
    
    // get the multiphasic material
    FEMultiphasic* pmb = m_pMat;
    const int nsol = (int)pmb->Solutes();
    vector< vector<double> > ct(nsol, vector<double>(FEElement::MAX_NODES));
    vector<int> sid(nsol);
    for (j=0; j<nsol; ++j) sid[j] = pmb->GetSolute(j)->GetSoluteDOF();
    
    // get the solid element
    FESolidElement& el = m_Elem[iel];
    
    // get the number of integration points
    nint = el.GaussPoints();
    
    // get the number of nodes
    neln = el.Nodes();
    
    // get the integration weights
    gw = el.GaussWeights();
    
    // get the nodal data
    for (j=0; j<neln; ++j)
    {
        FENode& node = mesh.Node(el.m_node[j]);
        r0[j] = node.m_r0;
        rt[j] = node.m_rt;
        if (el.m_bitfc.size()>0 && el.m_bitfc[j]) {
            pn[j] = (node.m_ID[m_dofQ] != -1) ? node.get(m_dofQ) : node.get(m_dofP);
            for (k=0; k<nsol; ++k)
                ct[k][j] = (node.m_ID[m_dofD + sid[k]] != -1) ? node.get(m_dofD + sid[k]) : node.get(m_dofC + sid[k]);
        }
        else {
            pn[j] = node.get(m_dofP);
            for (k=0; k<nsol; ++k)
                ct[k][j] = node.get(m_dofC + sid[k]);
        }
    }
    
    // loop over the integration points and calculate
    // the stress at the integration point
    for (n=0; n<nint; ++n)
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

        // multiphasic material point data
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        
        // update SBM referential densities
        pmb->UpdateSolidBoundMolecules(mp);
        
        // evaluate referential solid volume fraction
        ppt.m_phi0t = pmb->SolidReferentialVolumeFraction(mp);
        if (m_breset) ppt.m_phi0 = ppt.m_phi0t;

        // evaluate fluid pressure at gauss-point
        ppt.m_p = el.Evaluate(pn, n);
        
        // calculate the gradient of p at gauss-point
        ppt.m_gradp = gradient(el, pn, n);
        
        for (k=0; k<nsol; ++k) {
            // evaluate effective solute concentrations at gauss-point
            spt.m_c[k] = el.Evaluate(&ct[k][0], n);
            // calculate the gradient of c at gauss-point
            spt.m_gradc[k] = gradient(el, &ct[k][0], n);
        }
        
        // update the fluid and solute fluxes
        // and evaluate the actual fluid pressure and solute concentration
        ppt.m_w = pmb->FluidFlux(mp);
        spt.m_psi = pmb->ElectricPotential(mp);
        for (k=0; k<nsol; ++k) {
            spt.m_ca[k] = pmb->Concentration(mp,k);
            spt.m_j[k] = pmb->SoluteFlux(mp,k);
        }
        spt.m_cF = pmb->FixedChargeDensity(mp);
        ppt.m_pa = pmb->Pressure(mp);
        spt.m_Ie = pmb->CurrentDensity(mp);
        pmb->PartitionCoefficientFunctions(mp, spt.m_k, spt.m_dkdJ, spt.m_dkdc,
                                           spt.m_dkdr, spt.m_dkdJr, spt.m_dkdrc);

        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, GetFEModel()->GetTime());
        
        // calculate the solid stress at this material point
        ppt.m_ss = pmb->GetElasticMaterial()->Stress(mp);
        
        // evaluate the stress
        pt.m_s = pmb->Stress(mp);
        
        // evaluate the referential solid density
        spt.m_rhor = pmb->SolidReferentialApparentDensity(mp);
        
        // update chemical reaction element data
        for (int j=0; j<m_pMat->Reactions(); ++j)
            pmb->GetReaction(j)->UpdateElementData(mp);
        
    }
    if (m_breset) m_breset = false;
}
