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
#include "FEMultiphasicShellDomain.h"
#include "FEMultiphasicMultigeneration.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FEBioMech/FEBioMech.h>
#include <FECore/FELinearSystem.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
FEMultiphasicShellDomain::FEMultiphasicShellDomain(FEModel* pfem) : FESSIShellDomain(pfem), FEMultiphasicDomain(pfem), m_dofSU(pfem), m_dofR(pfem), m_dof(pfem)
{
    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofSU.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
        m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
    }
}

//-----------------------------------------------------------------------------
//! get the material (overridden from FEDomain)
FEMaterial* FEMultiphasicShellDomain::GetMaterial()
{ 
	return m_pMat; 
}

//-----------------------------------------------------------------------------
//! get the total dof
const FEDofList& FEMultiphasicShellDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
    m_pMat = dynamic_cast<FEMultiphasic*>(pmat);
    assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data.
void FEMultiphasicShellDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
    // get nodal DOFS
    const int nsol = m_pMat->Solutes();
    
    int N = el.Nodes();
    int ndpn = 8+2*nsol;
    lm.resize(N*ndpn);
    
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = m_pMesh->Node(n);
        
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[ndpn*i  ] = id[m_dofU[0]];
        lm[ndpn*i+1] = id[m_dofU[1]];
        lm[ndpn*i+2] = id[m_dofU[2]];
        
        // next the shell dofs
        lm[ndpn*i+3] = id[m_dofSU[0]];
        lm[ndpn*i+4] = id[m_dofSU[1]];
        lm[ndpn*i+5] = id[m_dofSU[2]];
        
        // now the pressure dofs
        lm[ndpn*i+6] = id[m_dofP];
        lm[ndpn*i+7] = id[m_dofQ];
        
        // concentration dofs in shell domain
        for (int k=0; k<nsol; ++k) {
            lm[ndpn*i+8+2*k] = id[m_dofC+m_pMat->GetSolute(k)->GetSoluteDOF()];
            lm[ndpn*i+9+2*k] = id[m_dofD+m_pMat->GetSolute(k)->GetSoluteDOF()];
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::UnpackMembraneLM(FEShellElement& el, vector<int>& lm)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // get the first material point
    FEMaterialPoint& mp = *el.GetMaterialPoint(0);
    FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());

    int nse = (int)ps.m_ide.size();
    int nsi = (int)ps.m_idi.size();
    vector<int> ide = ps.m_ide;
    vector<int> idi = ps.m_idi;

    int neln = el.Nodes();
    int ndpn = 8 + nse + nsi;
    lm.resize(neln*ndpn);

    for (int i=0; i<neln; ++i)
    {
        int n = el.m_node[i];
        
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[ndpn*i  ] = id[m_dofU[0]];
        lm[ndpn*i+1] = id[m_dofU[1]];
        lm[ndpn*i+2] = id[m_dofU[2]];
        
        // next the shell dofs
        lm[ndpn*i+3] = id[m_dofSU[0]];
        lm[ndpn*i+4] = id[m_dofSU[1]];
        lm[ndpn*i+5] = id[m_dofSU[2]];
        
        // now the pressure dofs
        lm[ndpn*i+6] = id[m_dofP];
        lm[ndpn*i+7] = id[m_dofQ];
        
        // concentration dofs
        for (int k=0; k<nse; ++k)
            lm[ndpn*i+8+k] = id[m_dofC + ide[k]];
        for (int k=0; k<nsi; ++k)
            lm[ndpn*i+8+nse+k] = id[m_dofD + idi[k]];
    }
}

//-----------------------------------------------------------------------------
//! build connectivity for matrix profile
void FEMultiphasicShellDomain::BuildMatrixProfile(FEGlobalMatrix& M)
{
    vector<int> elm;
    const int NE = Elements();
    for (int j = 0; j<NE; ++j)
    {
        FEElement& el = ElementRef(j);
        UnpackLM(el, elm);
        M.build_add(elm);
    }

    // if we have membrane reactions, build matrix profile for those
    if (m_pMat->MembraneReactions()) {
        for (int j = 0; j<NE; ++j)
        {
            FEShellElement& el = dynamic_cast<FEShellElement&>(ElementRef(j));
            UnpackMembraneLM(el, elm);
            M.build_add(elm);
        }
    }
}

//-----------------------------------------------------------------------------
bool FEMultiphasicShellDomain::Init()
{
    // initialize base class
	if (FESSIShellDomain::Init() == false) return false;
    
    // extract the initial concentrations of the solid-bound molecules
    const int nsbm = m_pMat->SBMs();
    const int nsol = m_pMat->Solutes();
    vector<double> sbmr(nsbm, 0);
    
    for (int i = 0; i<(int)m_Elem.size(); ++i)
    {
        UpdateShellMPData(i);
        
        // get the solid element
        FEShellElement& el = m_Elem[i];
        
        // get the number of integration points
        int nint = el.GaussPoints();
        
        // loop over the integration points
        for (int n = 0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            double h = el.Evaluate(el.m_ht, n);   // shell thickness

            // extract the initial apparent densities of the solid-bound molecules
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
            ps.m_sbmrp.assign(nsbm, 0);
            ps.m_sbmrhatp.assign(nsbm, 0);
            ps.m_sbmrhat.assign(nsbm,0);
            
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
            // check if this mixture includes membrane reactions
            int mreact = (int)m_pMat->MembraneReactions();
            if (mreact) {
                // for membrane reactions involving solid-bound molecules,
                // update their concentration
                // multiphasic material point data
                FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
                
                for (int isbm=0; isbm<nsbm; ++isbm) {
                    // combine the molar supplies from all the reactions
                    for (int k=0; k<mreact; ++k) {
                        double zetahat = m_pMat->GetMembraneReaction(k)->ReactionSupply(mp);
                        double v = m_pMat->GetMembraneReaction(k)->m_v[nsol+isbm];
                        // remember to convert from molar supply to referential mass supply
                        ps.m_sbmrhat[isbm] += pt.m_J/h*m_pMat->SBMMolarMass(isbm)*v*zetahat;
                    }
                }
            }
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::Activate()
{
    const int nsol = m_pMat->Solutes();
    
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
            for (int l=0; l<nsol; ++l)
            {
                int dofc = m_dofC + m_pMat->GetSolute(l)->GetSoluteDOF();
                node.set_active(dofc);
            }
            
            if (node.HasFlags(FENode::SHELL)) {
                node.set_active(m_dofQ);
                for (int l=0; l<nsol; ++l)
                {
                    int dofd = m_dofD + m_pMat->GetSolute(l)->GetSoluteDOF();
                    node.set_active(dofd);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::InitMaterialPoints()
{
    const int nsol = m_pMat->Solutes();
    const int nsbm = m_pMat->SBMs();
    FEMesh& m = *GetMesh();
    
    const int NE = FEElement::MAX_NODES;
    double p0[NE], q0[NE];
    vector< vector<double> > c0(nsol, vector<double>(NE));
    vector< vector<double> > d0(nsol, vector<double>(NE));
    vector<int> sid(nsol);
    for (int j = 0; j<nsol; ++j) sid[j] = m_pMat->GetSolute(j)->GetSoluteDOF();
    
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
            for (int isol = 0; isol<nsol; ++isol) {
                c0[isol][i] = m.Node(el.m_node[i]).get(m_dofC + sid[isol]);
                d0[isol][i] = m.Node(el.m_node[i]).get(m_dofD + sid[isol]);
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
            pt.m_p = evaluate(el, p0, q0, n);
            pt.m_gradp = gradient(el, p0, q0, n);
            pt.m_w = m_pMat->FluidFlux(mp);
            
            // initialize multiphasic solutes
            ps.m_nsol = nsol;
            ps.m_nsbm = nsbm;
            
            // initialize effective solute concentrations
            for (int isol = 0; isol<nsol; ++isol) {
                ps.m_c[isol] = evaluate(el, c0[isol], d0[isol], n);
                ps.m_gradc[isol] = gradient(el, c0[isol], d0[isol], n);
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
void FEMultiphasicShellDomain::Reset()
{
    // reset base class
    FESSIShellDomain::Reset();
    
    const int nsol = m_pMat->Solutes();
    const int nsbm = m_pMat->SBMs();
    vector<double> sbmr(nsbm,0);
    
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
            FEBiphasicMaterialPoint& pb = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());

            // extract the initial apparent densities of the solid-bound molecules
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
            ps.m_sbmrp.assign(nsbm, 0);
            ps.m_sbmrhatp.assign(nsbm, 0);
            ps.m_sbmrhat.assign(nsbm,0);
            
            // initialize referential solid volume fraction
            pb.m_phi0 = pb.m_phi0t = m_pMat->SolidReferentialVolumeFraction(mp);

            // reset chemical reaction element data
            ps.m_cri.clear();
            ps.m_crd.clear();
            for (int j=0; j<m_pMat->Reactions(); ++j)
                m_pMat->GetReaction(j)->ResetElementData(mp);
            
            // reset membrane reaction data
            if (m_pMat->MembraneReactions()) {
                ps.m_strain = 0;
                ps.m_pe = ps.m_pi = 0;
                ps.m_ide.clear();
                ps.m_idi.clear();
                ps.m_ce.clear();
                ps.m_ci.clear();
                for (int j=0; j<m_pMat->MembraneReactions(); ++j)
                    m_pMat->GetMembraneReaction(j)->ResetElementData(mp);
            }
        }
    }
    m_breset = true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    FESSIShellDomain::PreSolveUpdate(timeInfo);
    
    const int NE = FEElement::MAX_NODES;
    vec3d x0[NE], xt[NE], r0, rt;
    FEMesh& m = *GetMesh();
    
    for (size_t iel=0; iel<m_Elem.size(); ++iel)
    {
        FEShellElement& el = m_Elem[iel];
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
            for (int k=0; k<m_pMat->Solutes(); ++k) {
                ps.m_crp[k] = pe.m_J*m_pMat->Porosity(mp)*ps.m_ca[k];
            }
            
            // reset referential solid-bound molecule concentrations at previous time
            ps.m_sbmrp = ps.m_sbmr;
            ps.m_sbmrhatp = ps.m_sbmrhat;
            // Jansen initialization
            for (int k=0; k<ps.m_sbmrhat.size(); ++k)
                ps.m_sbmrhat[k] = - ps.m_sbmrhatp[k];
            
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
            
            // reset membrane reaction element data
            for (int j=0; j<m_pMat->MembraneReactions(); ++j)
                m_pMat->GetMembraneReaction(j)->InitializeElementData(mp);
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::InternalForces(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
    
    // get nodal DOFS
    int nsol = m_pMat->Solutes();
    int ndpn = 2*(4+nsol);
    
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = ndpn*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe, true);
    }
    
    MembraneReactionFluxes(R);
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEMultiphasicShellDomain::ElementInternalForce(FEShellElement& el, vector<double>& fe)
{
    int i, isol, n;
    
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
    
    const int nsol = m_pMat->Solutes();
    int ndpn = 2*(4+nsol);
    
    const int nreact = m_pMat->Reactions();
    const int mreact = m_pMat->MembraneReactions();
    
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
        
        // membrane reactions
        for (i=0; i<mreact; ++i) {
            FEMembraneReaction* pri = m_pMat->GetMembraneReaction(i);
            double zhat = pri->ReactionSupply(mp);
            phiwhat += phiw*pri->m_Vbar*zhat;
            for (isol=0; isol<nsol; ++isol)
                chat[isol] += phiw*zhat*pri->m_v[isol];
        }
        
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
            fe[ndpn*i  ] -= fu.x*detJt;
            fe[ndpn*i+1] -= fu.y*detJt;
            fe[ndpn*i+2] -= fu.z*detJt;
            fe[ndpn*i+3] -= fw.x*detJt;
            fe[ndpn*i+4] -= fw.y*detJt;
            fe[ndpn*i+5] -= fw.z*detJt;
            fe[ndpn*i+6] -= dt*(w*gradMu + (phiwhat - divv)*Mu)*detJt;
            fe[ndpn*i+7] -= dt*(w*gradMw + (phiwhat - divv)*Mw)*detJt;
            for (isol=0; isol<nsol; ++isol) {
                fe[ndpn*i+8+2*isol] -= dt*(gradMu*(j[isol]+je*m_pMat->m_penalty)
                                           + Mu*(chat[isol] - (phiw*spt.m_ca[isol] - spt.m_crp[isol]/J)/dt)
                                           )*detJt;
                fe[ndpn*i+9+2*isol] -= dt*(gradMw*(j[isol]+je*m_pMat->m_penalty)
                                           + Mw*(chat[isol] - (phiw*spt.m_ca[isol] - spt.m_crp[isol]/J)/dt)
                                           )*detJt;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::InternalForcesSS(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
    
    // get nodal DOFS
    int nsol = m_pMat->Solutes();
    int ndpn = 2*(4+nsol);
    
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = ndpn*el.Nodes();
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

void FEMultiphasicShellDomain::ElementInternalForceSS(FEShellElement& el, vector<double>& fe)
{
    int i, isol, n;
    
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
    
    vec3d gcnt[3];
    
    const int nsol = m_pMat->Solutes();
    int ndpn = 2*(4+nsol);
    
    const int nreact = m_pMat->Reactions();
    const int mreact = m_pMat->MembraneReactions();

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
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
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
        
        // membrane reactions
        for (i=0; i<mreact; ++i) {
            FEMembraneReaction* pri = m_pMat->GetMembraneReaction(i);
            double zhat = pri->ReactionSupply(mp);
            phiwhat += phiw*pri->m_Vbar*zhat;
            for (isol=0; isol<nsol; ++isol)
                chat[isol] += phiw*zhat*pri->m_v[isol];
        }
        
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
            fe[ndpn*i  ] -= fu.x*detJt;
            fe[ndpn*i+1] -= fu.y*detJt;
            fe[ndpn*i+2] -= fu.z*detJt;
            fe[ndpn*i+3] -= fw.x*detJt;
            fe[ndpn*i+4] -= fw.y*detJt;
            fe[ndpn*i+5] -= fw.z*detJt;
            fe[ndpn*i+6] -= dt*(w*gradMu + Mu*phiwhat)*detJt;
            fe[ndpn*i+7] -= dt*(w*gradMw + Mw*phiwhat)*detJt;
            for (isol=0; isol<nsol; ++isol) {
                fe[ndpn*i+8+2*isol] -= dt*(gradMu*(j[isol]+je*m_pMat->m_penalty)
                                         + Mu*phiw*chat[isol]
                                         )*detJt;
                fe[ndpn*i+9+2*isol] -= dt*(gradMw*(j[isol]+je*m_pMat->m_penalty)
                                           + Mw*phiw*chat[isol]
                                           )*detJt;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::StiffnessMatrix(FELinearSystem& LS, bool bsymm)
{
    const int nsol = m_pMat->Solutes();
    int ndpn = 2*(4+nsol);
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FEShellElement& el = m_Elem[iel];

        // element stiffness matrix
		FEElementMatrix ke(el);
		int neln = el.Nodes();
        int ndof = neln*ndpn;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementMultiphasicStiffness(el, ke, bsymm);

		// get lm vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
    
    MembraneReactionStiffnessMatrix(LS);
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::StiffnessMatrixSS(FELinearSystem& LS, bool bsymm)
{
    const int nsol = m_pMat->Solutes();
    int ndpn = 2*(4+nsol);
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FEShellElement& el = m_Elem[iel];

        // element stiffness matrix
		FEElementMatrix ke(el);
        int neln = el.Nodes();
        int ndof = neln*ndpn;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementMultiphasicStiffnessSS(el, ke, bsymm);

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
bool FEMultiphasicShellDomain::ElementMultiphasicStiffness(FEShellElement& el, matrix& ke, bool bsymm)
{
    int i, j, isol, jsol, n, ireact, isbm;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
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
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    const int nsol = m_pMat->Solutes();
    int ndpn = 2*(4+nsol);
    
    const int nsbm   = m_pMat->SBMs();
    const int nreact = m_pMat->Reactions();
    const int mreact = m_pMat->MembraneReactions();

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
        for (isol=0; isol<nsol; ++isol)
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
        for (i=0; i<nreact; ++i)
            Phie += m_pMat->GetReaction(i)->m_Vbar*(I*m_pMat->GetReaction(i)->ReactionSupply(mp)
                                                    +m_pMat->GetReaction(i)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
        
        // membrane reactions
        for (i=0; i<mreact; ++i)
            Phie += m_pMat->GetMembraneReaction(i)->m_Vbar*mat3dd(m_pMat->GetMembraneReaction(i)->ReactionSupply(mp)
                                                    +m_pMat->GetMembraneReaction(i)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
        
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
            
            // chemical reactions
            dchatde[isol].zero();
            for (ireact=0; ireact<nreact; ++ireact) {
                dchatde[isol] += m_pMat->GetReaction(ireact)->m_v[isol]
                *(I*m_pMat->GetReaction(ireact)->ReactionSupply(mp)
                  +m_pMat->GetReaction(ireact)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
                Phic[isol] += phiw*m_pMat->GetReaction(ireact)->m_Vbar
                *m_pMat->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp, isol);
            }
            
            // membrane reactions
            for (ireact=0; ireact<mreact; ++ireact) {
                dchatde[isol] += m_pMat->GetMembraneReaction(ireact)->m_v[isol]
                *mat3dd(m_pMat->GetMembraneReaction(ireact)->ReactionSupply(mp)
                  +m_pMat->GetMembraneReaction(ireact)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
                Phic[isol] += phiw*m_pMat->GetMembraneReaction(ireact)->m_Vbar
                *m_pMat->GetMembraneReaction(ireact)->Tangent_ReactionSupply_Concentration(mp, isol);
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
        vec3d vtmp,gp,qpu, qpw;
        vector<vec3d> gc(nsol),qcu(nsol),qcw(nsol),wc(nsol),wd(nsol),jce(nsol),jde(nsol);
        vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
        vector< vector<vec3d> > jd(nsol, vector<vec3d>(nsol));
        mat3d wu, ww, jue, jwe;
        vector<mat3d> ju(nsol), jw(nsol);
        vector< vector<double> > qcc(nsol, vector<double>(nsol));
        vector< vector<double> > qcd(nsol, vector<double>(nsol));
        vector< vector<double> > dchatdc(nsol, vector<double>(nsol));
        double sum;
        mat3ds De;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradMu[i]*(s*gradMu[j])) + vdotTdotv(gradMu[i], C, gradMu[j]))*detJ;
                mat3d Kuw = (mat3dd(gradMu[i]*(s*gradMw[j])) + vdotTdotv(gradMu[i], C, gradMw[j]))*detJ;
                mat3d Kwu = (mat3dd(gradMw[i]*(s*gradMu[j])) + vdotTdotv(gradMw[i], C, gradMu[j]))*detJ;
                mat3d Kww = (mat3dd(gradMw[i]*(s*gradMw[j])) + vdotTdotv(gradMw[i], C, gradMw[j]))*detJ;

                ke[ndpn*i  ][ndpn*j  ] += Kuu[0][0]; ke[ndpn*i  ][ndpn*j+1] += Kuu[0][1]; ke[ndpn*i  ][ndpn*j+2] += Kuu[0][2];
                ke[ndpn*i+1][ndpn*j  ] += Kuu[1][0]; ke[ndpn*i+1][ndpn*j+1] += Kuu[1][1]; ke[ndpn*i+1][ndpn*j+2] += Kuu[1][2];
                ke[ndpn*i+2][ndpn*j  ] += Kuu[2][0]; ke[ndpn*i+2][ndpn*j+1] += Kuu[2][1]; ke[ndpn*i+2][ndpn*j+2] += Kuu[2][2];
                
                ke[ndpn*i  ][ndpn*j+3] += Kuw[0][0]; ke[ndpn*i  ][ndpn*j+4] += Kuw[0][1]; ke[ndpn*i  ][ndpn*j+5] += Kuw[0][2];
                ke[ndpn*i+1][ndpn*j+3] += Kuw[1][0]; ke[ndpn*i+1][ndpn*j+4] += Kuw[1][1]; ke[ndpn*i+1][ndpn*j+5] += Kuw[1][2];
                ke[ndpn*i+2][ndpn*j+3] += Kuw[2][0]; ke[ndpn*i+2][ndpn*j+4] += Kuw[2][1]; ke[ndpn*i+2][ndpn*j+5] += Kuw[2][2];
                
                ke[ndpn*i+3][ndpn*j  ] += Kwu[0][0]; ke[ndpn*i+3][ndpn*j+1] += Kwu[0][1]; ke[ndpn*i+3][ndpn*j+2] += Kwu[0][2];
                ke[ndpn*i+4][ndpn*j  ] += Kwu[1][0]; ke[ndpn*i+4][ndpn*j+1] += Kwu[1][1]; ke[ndpn*i+4][ndpn*j+2] += Kwu[1][2];
                ke[ndpn*i+5][ndpn*j  ] += Kwu[2][0]; ke[ndpn*i+5][ndpn*j+1] += Kwu[2][1]; ke[ndpn*i+5][ndpn*j+2] += Kwu[2][2];
                
                ke[ndpn*i+3][ndpn*j+3] += Kww[0][0]; ke[ndpn*i+3][ndpn*j+4] += Kww[0][1]; ke[ndpn*i+3][ndpn*j+5] += Kww[0][2];
                ke[ndpn*i+4][ndpn*j+3] += Kww[1][0]; ke[ndpn*i+4][ndpn*j+4] += Kww[1][1]; ke[ndpn*i+4][ndpn*j+5] += Kww[1][2];
                ke[ndpn*i+5][ndpn*j+3] += Kww[2][0]; ke[ndpn*i+5][ndpn*j+4] += Kww[2][1]; ke[ndpn*i+5][ndpn*j+5] += Kww[2][2];
                
                // calculate the kpu matrix
                gp = vec3d(0,0,0);
                for (isol=0; isol<nsol; ++isol) gp += (D[isol]*gradc[isol])*(kappa[isol]/D0[isol]);
                gp = gradp+gp*(R*T);
                wu = vdotTdotv(-gp, dKedE, gradMu[j]);
                ww = vdotTdotv(-gp, dKedE, gradMw[j]);
                for (isol=0; isol<nsol; ++isol) {
                    wu += (((Ke*(D[isol]*gradc[isol])) & gradMu[j])*(J*dkdJ[isol] - kappa[isol])
                           +Ke*(2*kappa[isol]*(gradMu[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
                    + (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradMu[j]))*(-kappa[isol]*R*T/D0[isol]);
                    ww += (((Ke*(D[isol]*gradc[isol])) & gradMw[j])*(J*dkdJ[isol] - kappa[isol])
                           +Ke*(2*kappa[isol]*(gradMw[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
                    + (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradMw[j]))*(-kappa[isol]*R*T/D0[isol]);
                }
                qpu = -gradMu[j]*(1.0/dt);
                qpw = -gradMw[j]*(1.0/dt);
                vec3d kpu = (wu.transpose()*gradMu[i] + (qpu + Phie*gradMu[j])*Mu[i])*(detJ*dt);
                vec3d kpw = (ww.transpose()*gradMu[i] + (qpw + Phie*gradMw[j])*Mu[i])*(detJ*dt);
                vec3d kqu = (wu.transpose()*gradMw[i] + (qpu + Phie*gradMu[j])*Mw[i])*(detJ*dt);
                vec3d kqw = (ww.transpose()*gradMw[i] + (qpw + Phie*gradMw[j])*Mw[i])*(detJ*dt);
                ke[ndpn*i+6][ndpn*j  ] += kpu.x; ke[ndpn*i+6][ndpn*j+1] += kpu.y; ke[ndpn*i+6][ndpn*j+2] += kpu.z;
                ke[ndpn*i+6][ndpn*j+3] += kpw.x; ke[ndpn*i+6][ndpn*j+4] += kpw.y; ke[ndpn*i+6][ndpn*j+5] += kpw.z;
                ke[ndpn*i+7][ndpn*j  ] += kqu.x; ke[ndpn*i+7][ndpn*j+1] += kqu.y; ke[ndpn*i+7][ndpn*j+2] += kqu.z;
                ke[ndpn*i+7][ndpn*j+3] += kqw.x; ke[ndpn*i+7][ndpn*j+4] += kqw.y; ke[ndpn*i+7][ndpn*j+5] += kqw.z;
                
                // calculate the kup matrix
                vec3d kup = gradMu[i]*(-Mu[j]*detJ);
                vec3d kuq = gradMu[i]*(-Mw[j]*detJ);
                vec3d kwp = gradMw[i]*(-Mu[j]*detJ);
                vec3d kwq = gradMw[i]*(-Mw[j]*detJ);

                ke[ndpn*i  ][ndpn*j+6] += kup.x; ke[ndpn*i  ][ndpn*j+7] += kuq.x;
                ke[ndpn*i+1][ndpn*j+6] += kup.y; ke[ndpn*i+1][ndpn*j+7] += kuq.y;
                ke[ndpn*i+2][ndpn*j+6] += kup.z; ke[ndpn*i+2][ndpn*j+7] += kuq.z;
                
                ke[ndpn*i+3][ndpn*j+6] += kwp.x; ke[ndpn*i+3][ndpn*j+7] += kwq.x;
                ke[ndpn*i+4][ndpn*j+6] += kwp.y; ke[ndpn*i+4][ndpn*j+7] += kwq.y;
                ke[ndpn*i+5][ndpn*j+6] += kwp.z; ke[ndpn*i+5][ndpn*j+7] += kwq.z;
                
                // calculate the kpp matrix
                ke[ndpn*i+6][ndpn*j+6] += (Mu[i]*Mu[j]*Phip - gradMu[i]*(Ke*gradMu[j]))*(detJ*dt);
                ke[ndpn*i+6][ndpn*j+7] += (Mu[i]*Mw[j]*Phip - gradMu[i]*(Ke*gradMw[j]))*(detJ*dt);
                ke[ndpn*i+7][ndpn*j+6] += (Mw[i]*Mu[j]*Phip - gradMw[i]*(Ke*gradMu[j]))*(detJ*dt);
                ke[ndpn*i+7][ndpn*j+7] += (Mw[i]*Mw[j]*Phip - gradMw[i]*(Ke*gradMw[j]))*(detJ*dt);
                
                // calculate kcu matrix data
                jue.zero(); jwe.zero();
                De.zero();
                for (isol=0; isol<nsol; ++isol) {
                    gc[isol] = -gradc[isol]*phiw + w*c[isol]/D0[isol];
                    ju[isol] = ((D[isol]*gc[isol]) & gradMu[j])*(J*dkdJ[isol])
                    + vdotTdotv(gc[isol], dDdE[isol], gradMu[j])*kappa[isol]
                    + (((D[isol]*gradc[isol]) & gradMu[j])*(-phis)
                       +(D[isol]*((gradMu[j]*w)*2) - ((D[isol]*w) & gradMu[j]))*c[isol]/D0[isol]
                       )*kappa[isol]
                    +D[isol]*wu*(kappa[isol]*c[isol]/D0[isol]);
                    jw[isol] = ((D[isol]*gc[isol]) & gradMw[j])*(J*dkdJ[isol])
                    + vdotTdotv(gc[isol], dDdE[isol], gradMw[j])*kappa[isol]
                    + (((D[isol]*gradc[isol]) & gradMw[j])*(-phis)
                       +(D[isol]*((gradMw[j]*w)*2) - ((D[isol]*w) & gradMw[j]))*c[isol]/D0[isol]
                       )*kappa[isol]
                    +D[isol]*ww*(kappa[isol]*c[isol]/D0[isol]);
                    jue += ju[isol]*z[isol];
                    jwe += jw[isol]*z[isol];
                    De += D[isol]*(z[isol]*kappa[isol]*c[isol]/D0[isol]);
                    qcu[isol] = qpu*(c[isol]*(kappa[isol]+J*phiw*dkdJ[isol]));
                    qcw[isol] = qpw*(c[isol]*(kappa[isol]+J*phiw*dkdJ[isol]));
                    
                    // chemical reactions
                    for (ireact=0; ireact<nreact; ++ireact) {
                        double sum1 = 0;
                        double sum2 = 0;
                        for (isbm=0; isbm<nsbm; ++isbm) {
                            sum1 += m_pMat->SBMMolarMass(isbm)*m_pMat->GetReaction(ireact)->m_v[nsol+isbm]*
                            ((J-phi0)*dkdr[isol][isbm]-kappa[isol]/m_pMat->SBMDensity(isbm));
                            sum2 += m_pMat->SBMMolarMass(isbm)*m_pMat->GetReaction(ireact)->m_v[nsol+isbm]*
                            (dkdr[isol][isbm]+(J-phi0)*dkdJr[isol][isbm]-dkdJ[isol]/m_pMat->SBMDensity(isbm));
                        }
                        double zhat = m_pMat->GetReaction(ireact)->ReactionSupply(mp);
                        mat3dd zhatI(zhat);
                        mat3ds dzde = m_pMat->GetReaction(ireact)->Tangent_ReactionSupply_Strain(mp);
                        qcu[isol] -= ((zhatI+dzde*(J-phi0))*gradMu[j])*(sum1*c[isol])
                        +gradMu[j]*(c[isol]*(J-phi0)*sum2*zhat);
                        qcw[isol] -= ((zhatI+dzde*(J-phi0))*gradMw[j])*(sum1*c[isol])
                        +gradMw[j]*(c[isol]*(J-phi0)*sum2*zhat);
                    }
                    
                    // membrane reactions
                    for (ireact=0; ireact<mreact; ++ireact) {
                        double sum1 = 0;
                        double sum2 = 0;
                        for (isbm=0; isbm<nsbm; ++isbm) {
                            sum1 += m_pMat->SBMMolarMass(isbm)*m_pMat->GetMembraneReaction(ireact)->m_v[nsol+isbm]*
                            ((J-phi0)*dkdr[isol][isbm]-kappa[isol]/m_pMat->SBMDensity(isbm));
                            sum2 += m_pMat->SBMMolarMass(isbm)*m_pMat->GetMembraneReaction(ireact)->m_v[nsol+isbm]*
                            (dkdr[isol][isbm]+(J-phi0)*dkdJr[isol][isbm]-dkdJ[isol]/m_pMat->SBMDensity(isbm));
                        }
                        double zhat = m_pMat->GetMembraneReaction(ireact)->ReactionSupply(mp);
                        mat3dd zhatI(zhat);
                        mat3ds dzde = mat3dd(m_pMat->GetMembraneReaction(ireact)->Tangent_ReactionSupply_Strain(mp));
                        qcu[isol] -= ((zhatI+dzde*(J-phi0))*gradMu[j])*(sum1*c[isol])
                        +gradMu[j]*(c[isol]*(J-phi0)*sum2*zhat);
                        qcw[isol] -= ((zhatI+dzde*(J-phi0))*gradMw[j])*(sum1*c[isol])
                        +gradMw[j]*(c[isol]*(J-phi0)*sum2*zhat);
                    }
                }
                
                for (isol=0; isol<nsol; ++isol) {
                    
                    // calculate the kcu matrix
                    vec3d kcu = ((ju[isol]+jue*penalty).transpose()*gradMu[i]
                            + (qcu[isol] + dchatde[isol]*gradMu[j])*Mu[i])*(detJ*dt);
                    vec3d kcw = ((jw[isol]+jwe*penalty).transpose()*gradMu[i]
                                 + (qcw[isol] + dchatde[isol]*gradMw[j])*Mu[i])*(detJ*dt);
                    vec3d kdu = ((ju[isol]+jue*penalty).transpose()*gradMw[i]
                                 + (qcu[isol] + dchatde[isol]*gradMu[j])*Mw[i])*(detJ*dt);
                    vec3d kdw = ((jw[isol]+jwe*penalty).transpose()*gradMw[i]
                                 + (qcw[isol] + dchatde[isol]*gradMw[j])*Mw[i])*(detJ*dt);
                    ke[ndpn*i+8+2*isol][ndpn*j  ] += kcu.x; ke[ndpn*i+8+2*isol][ndpn*j+1] += kcu.y; ke[ndpn*i+8+2*isol][ndpn*j+2] += kcu.z;
                    ke[ndpn*i+8+2*isol][ndpn*j+3] += kcw.x; ke[ndpn*i+8+2*isol][ndpn*j+4] += kcw.y; ke[ndpn*i+8+2*isol][ndpn*j+5] += kcw.z;
                    ke[ndpn*i+9+2*isol][ndpn*j  ] += kdu.x; ke[ndpn*i+9+2*isol][ndpn*j+1] += kdu.y; ke[ndpn*i+9+2*isol][ndpn*j+2] += kdu.z;
                    ke[ndpn*i+9+2*isol][ndpn*j+3] += kdw.x; ke[ndpn*i+9+2*isol][ndpn*j+4] += kdw.y; ke[ndpn*i+9+2*isol][ndpn*j+5] += kdw.z;
                    
                    // calculate the kcp matrix
                    ke[ndpn*i+8+2*isol][ndpn*j+6] -= (gradMu[i]*(
                                                                 (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                                  +De*penalty)
                                                                 *(Ke*gradMu[j])
                                                                 ))*(detJ*dt);
                    ke[ndpn*i+8+2*isol][ndpn*j+7] -= (gradMu[i]*(
                                                                 (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                                  +De*penalty)
                                                                 *(Ke*gradMw[j])
                                                                 ))*(detJ*dt);
                    ke[ndpn*i+9+2*isol][ndpn*j+6] -= (gradMw[i]*(
                                                                 (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                                  +De*penalty)
                                                                 *(Ke*gradMu[j])
                                                                 ))*(detJ*dt);
                    ke[ndpn*i+9+2*isol][ndpn*j+7] -= (gradMw[i]*(
                                                                 (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                                  +De*penalty)
                                                                 *(Ke*gradMw[j])
                                                                 ))*(detJ*dt);
                    
                    // calculate the kuc matrix
                    sum = 0;
                    for (jsol=0; jsol<nsol; ++jsol)
                        sum += c[jsol]*(dodc[isol]*kappa[jsol]+osmc*dkdc[jsol][isol]);
                    vec3d kuc = (dTdc[isol]*gradMu[i] - gradMu[i]*(R*T*(osmc*kappa[isol]+sum)))*Mu[j]*detJ;
                    vec3d kud = (dTdc[isol]*gradMu[i] - gradMu[i]*(R*T*(osmc*kappa[isol]+sum)))*Mw[j]*detJ;
                    vec3d kwc = (dTdc[isol]*gradMw[i] - gradMw[i]*(R*T*(osmc*kappa[isol]+sum)))*Mu[j]*detJ;
                    vec3d kwd = (dTdc[isol]*gradMw[i] - gradMw[i]*(R*T*(osmc*kappa[isol]+sum)))*Mw[j]*detJ;

                    ke[ndpn*i  ][ndpn*j+8+2*isol] += kuc.x; ke[ndpn*i  ][ndpn*j+9+2*isol] += kud.x;
                    ke[ndpn*i+1][ndpn*j+8+2*isol] += kuc.y; ke[ndpn*i+1][ndpn*j+9+2*isol] += kud.y;
                    ke[ndpn*i+2][ndpn*j+8+2*isol] += kuc.z; ke[ndpn*i+2][ndpn*j+9+2*isol] += kud.z;
                    
                    ke[ndpn*i+3][ndpn*j+8+2*isol] += kwc.x; ke[ndpn*i+3][ndpn*j+9+2*isol] += kwd.x;
                    ke[ndpn*i+4][ndpn*j+8+2*isol] += kwc.y; ke[ndpn*i+4][ndpn*j+9+2*isol] += kwd.y;
                    ke[ndpn*i+5][ndpn*j+8+2*isol] += kwc.z; ke[ndpn*i+5][ndpn*j+9+2*isol] += kwd.z;
                    
                    // calculate the kpc matrix
                    vtmp = vec3d(0,0,0);
                    for (jsol=0; jsol<nsol; ++jsol)
                        vtmp += (D[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol])
                                 +dDdc[jsol][isol]*kappa[jsol])/D0[jsol]*gradc[jsol];
                    wc[isol] = (dKedc[isol]*gp)*(-Mu[j])
                    -Ke*((D[isol]*gradMu[j])*(kappa[isol]/D0[isol])+vtmp*Mu[j])*(R*T);
                    wd[isol] = (dKedc[isol]*gp)*(-Mw[j])
                    -Ke*((D[isol]*gradMw[j])*(kappa[isol]/D0[isol])+vtmp*Mw[j])*(R*T);

                    ke[ndpn*i+6][ndpn*j+8+2*isol] += (gradMu[i]*wc[isol])*(detJ*dt);
                    ke[ndpn*i+6][ndpn*j+9+2*isol] += (gradMu[i]*wd[isol])*(detJ*dt);
                    ke[ndpn*i+7][ndpn*j+8+2*isol] += (gradMw[i]*wc[isol])*(detJ*dt);
                    ke[ndpn*i+7][ndpn*j+9+2*isol] += (gradMw[i]*wd[isol])*(detJ*dt);
                    
                }
                
                // calculate data for the kcc matrix
                jce.assign(nsol, vec3d(0,0,0));
                jde.assign(nsol, vec3d(0,0,0));
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        if (jsol != isol) {
                            jc[isol][jsol] =
                            ((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*Mu[j]
                            +(D[isol]*(w*(-Mu[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            jd[isol][jsol] =
                            ((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*Mw[j]
                            +(D[isol]*(w*(-Mw[j]*dD0dc[isol][jsol]/D0[isol])+wd[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            
                            qcc[isol][jsol] = -Mu[j]*phiw/dt*c[isol]*dkdc[isol][jsol];
                            qcd[isol][jsol] = -Mw[j]*phiw/dt*c[isol]*dkdc[isol][jsol];
                        }
                        else {
                            jc[isol][jsol] = (D[isol]*(gradMu[j]*(-phiw)+w*(Mu[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*Mu[j]
                            +(D[isol]*(w*(-Mu[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            jd[isol][jsol] = (D[isol]*(gradMw[j]*(-phiw)+w*(Mw[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*Mw[j]
                            +(D[isol]*(w*(-Mw[j]*dD0dc[isol][jsol]/D0[isol])+wd[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            
                            qcc[isol][jsol] = -Mu[j]*phiw/dt*(c[isol]*dkdc[isol][jsol] + kappa[isol]);
                            qcd[isol][jsol] = -Mw[j]*phiw/dt*(c[isol]*dkdc[isol][jsol] + kappa[isol]);
                        }
                        jce[jsol] += jc[isol][jsol]*z[isol];
                        jde[jsol] += jd[isol][jsol]*z[isol];
                        
                        // chemical reactions
                        dchatdc[isol][jsol] = 0;
                        for (ireact=0; ireact<nreact; ++ireact) {
                            dchatdc[isol][jsol] += m_pMat->GetReaction(ireact)->m_v[isol]
                            *m_pMat->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp,jsol);
                            double sum1 = 0;
                            double sum2 = 0;
                            for (isbm=0; isbm<nsbm; ++isbm) {
                                sum1 += m_pMat->SBMMolarMass(isbm)*m_pMat->GetReaction(ireact)->m_v[nsol+isbm]*
                                ((J-phi0)*dkdr[isol][isbm]-kappa[isol]/m_pMat->SBMDensity(isbm));
                                sum2 += m_pMat->SBMMolarMass(isbm)*m_pMat->GetReaction(ireact)->m_v[nsol+isbm]*
                                ((J-phi0)*dkdrc[isol][isbm][jsol]-dkdc[isol][jsol]/m_pMat->SBMDensity(isbm));
                            }
                            double zhat = m_pMat->GetReaction(ireact)->ReactionSupply(mp);
                            double dzdc = m_pMat->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp, jsol);
                            if (jsol != isol) {
                                qcc[isol][jsol] -= Mu[j]*phiw*c[isol]*(dzdc*sum1+zhat*sum2);
                                qcd[isol][jsol] -= Mw[j]*phiw*c[isol]*(dzdc*sum1+zhat*sum2);
                            }
                            else {
                                qcc[isol][jsol] -= Mu[j]*phiw*((zhat+c[isol]*dzdc)*sum1+c[isol]*zhat*sum2);
                                qcd[isol][jsol] -= Mw[j]*phiw*((zhat+c[isol]*dzdc)*sum1+c[isol]*zhat*sum2);
                            }
                        }
                        
                        // membrane reactions
                        for (ireact=0; ireact<mreact; ++ireact) {
                            dchatdc[isol][jsol] += m_pMat->GetMembraneReaction(ireact)->m_v[isol]
                            *m_pMat->GetMembraneReaction(ireact)->Tangent_ReactionSupply_Concentration(mp,jsol);
                            double sum1 = 0;
                            double sum2 = 0;
                            for (isbm=0; isbm<nsbm; ++isbm) {
                                sum1 += m_pMat->SBMMolarMass(isbm)*m_pMat->GetMembraneReaction(ireact)->m_v[nsol+isbm]*
                                ((J-phi0)*dkdr[isol][isbm]-kappa[isol]/m_pMat->SBMDensity(isbm));
                                sum2 += m_pMat->SBMMolarMass(isbm)*m_pMat->GetMembraneReaction(ireact)->m_v[nsol+isbm]*
                                ((J-phi0)*dkdrc[isol][isbm][jsol]-dkdc[isol][jsol]/m_pMat->SBMDensity(isbm));
                            }
                            double zhat = m_pMat->GetMembraneReaction(ireact)->ReactionSupply(mp);
                            double dzdc = m_pMat->GetMembraneReaction(ireact)->Tangent_ReactionSupply_Concentration(mp, jsol);
                            if (jsol != isol) {
                                qcc[isol][jsol] -= Mu[j]*phiw*c[isol]*(dzdc*sum1+zhat*sum2);
                                qcd[isol][jsol] -= Mw[j]*phiw*c[isol]*(dzdc*sum1+zhat*sum2);
                            }
                            else {
                                qcc[isol][jsol] -= Mu[j]*phiw*((zhat+c[isol]*dzdc)*sum1+c[isol]*zhat*sum2);
                                qcd[isol][jsol] -= Mw[j]*phiw*((zhat+c[isol]*dzdc)*sum1+c[isol]*zhat*sum2);
                            }
                        }
                    }
                }
                
                // calculate the kcc matrix
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        ke[ndpn*i+8+2*isol][ndpn*j+8+2*jsol] += (gradMu[i]*(jc[isol][jsol]+jce[jsol]*penalty)
                                                             + Mu[i]*(qcc[isol][jsol]
                                                                     + Mu[j]*phiw*dchatdc[isol][jsol]))*(detJ*dt);
                        ke[ndpn*i+8+2*isol][ndpn*j+9+2*jsol] += (gradMu[i]*(jd[isol][jsol]+jde[jsol]*penalty)
                                                                 + Mu[i]*(qcd[isol][jsol]
                                                                          + Mw[j]*phiw*dchatdc[isol][jsol]))*(detJ*dt);
                        ke[ndpn*i+9+2*isol][ndpn*j+8+2*jsol] += (gradMw[i]*(jc[isol][jsol]+jce[jsol]*penalty)
                                                                 + Mw[i]*(qcc[isol][jsol]
                                                                          + Mu[j]*phiw*dchatdc[isol][jsol]))*(detJ*dt);
                        ke[ndpn*i+9+2*isol][ndpn*j+9+2*jsol] += (gradMw[i]*(jd[isol][jsol]+jde[jsol]*penalty)
                                                                 + Mw[i]*(qcd[isol][jsol]
                                                                          + Mw[j]*phiw*dchatdc[isol][jsol]))*(detJ*dt);
                    }
                }
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
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
//! calculates element stiffness matrix for element iel
//! for steady-state response (zero solid velocity, zero time derivative of
//! solute concentration)
//!
bool FEMultiphasicShellDomain::ElementMultiphasicStiffnessSS(FEShellElement& el, matrix& ke, bool bsymm)
{
    int i, j, isol, jsol, n, ireact;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
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
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    const int nsol = m_pMat->Solutes();
    int ndpn = 2*(4+nsol);
    
    const int nreact = m_pMat->Reactions();
    const int mreact = m_pMat->MembraneReactions();

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
        mat3ds Phie; Phie.zero();
        double Phip = 0;
        vector<double> Phic(nsol,0);
        if (m_pMat->GetSolventSupply()) {
            Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
            Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
        }
        
        // chemical reactions
        for (i=0; i<nreact; ++i)
            Phie += m_pMat->GetReaction(i)->m_Vbar*(I*m_pMat->GetReaction(i)->ReactionSupply(mp)
                                                    +m_pMat->GetReaction(i)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
        
        // membrane reactions
        for (i=0; i<mreact; ++i)
            Phie += m_pMat->GetReaction(i)->m_Vbar*mat3dd(m_pMat->GetMembraneReaction(i)->ReactionSupply(mp)
                                                    +m_pMat->GetMembraneReaction(i)->Tangent_ReactionSupply_Strain(mp)*(J*phiw));
        
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
        vec3d vtmp,gp,qpu, qpw;
        vector<vec3d> gc(nsol),wc(nsol),wd(nsol),jce(nsol),jde(nsol);
        vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
        vector< vector<vec3d> > jd(nsol, vector<vec3d>(nsol));
        mat3d wu, ww, jue, jwe;
        vector<mat3d> ju(nsol), jw(nsol);
        vector< vector<double> > dchatdc(nsol, vector<double>(nsol));
        double sum;
        mat3ds De;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradMu[i]*(s*gradMu[j])) + vdotTdotv(gradMu[i], C, gradMu[j]))*detJ;
                mat3d Kuw = (mat3dd(gradMu[i]*(s*gradMw[j])) + vdotTdotv(gradMu[i], C, gradMw[j]))*detJ;
                mat3d Kwu = (mat3dd(gradMw[i]*(s*gradMu[j])) + vdotTdotv(gradMw[i], C, gradMu[j]))*detJ;
                mat3d Kww = (mat3dd(gradMw[i]*(s*gradMw[j])) + vdotTdotv(gradMw[i], C, gradMw[j]))*detJ;
                
                ke[ndpn*i  ][ndpn*j  ] += Kuu[0][0]; ke[ndpn*i  ][ndpn*j+1] += Kuu[0][1]; ke[ndpn*i  ][ndpn*j+2] += Kuu[0][2];
                ke[ndpn*i+1][ndpn*j  ] += Kuu[1][0]; ke[ndpn*i+1][ndpn*j+1] += Kuu[1][1]; ke[ndpn*i+1][ndpn*j+2] += Kuu[1][2];
                ke[ndpn*i+2][ndpn*j  ] += Kuu[2][0]; ke[ndpn*i+2][ndpn*j+1] += Kuu[2][1]; ke[ndpn*i+2][ndpn*j+2] += Kuu[2][2];
                
                ke[ndpn*i  ][ndpn*j+3] += Kuw[0][0]; ke[ndpn*i  ][ndpn*j+4] += Kuw[0][1]; ke[ndpn*i  ][ndpn*j+5] += Kuw[0][2];
                ke[ndpn*i+1][ndpn*j+3] += Kuw[1][0]; ke[ndpn*i+1][ndpn*j+4] += Kuw[1][1]; ke[ndpn*i+1][ndpn*j+5] += Kuw[1][2];
                ke[ndpn*i+2][ndpn*j+3] += Kuw[2][0]; ke[ndpn*i+2][ndpn*j+4] += Kuw[2][1]; ke[ndpn*i+2][ndpn*j+5] += Kuw[2][2];
                
                ke[ndpn*i+3][ndpn*j  ] += Kwu[0][0]; ke[ndpn*i+3][ndpn*j+1] += Kwu[0][1]; ke[ndpn*i+3][ndpn*j+2] += Kwu[0][2];
                ke[ndpn*i+4][ndpn*j  ] += Kwu[1][0]; ke[ndpn*i+4][ndpn*j+1] += Kwu[1][1]; ke[ndpn*i+4][ndpn*j+2] += Kwu[1][2];
                ke[ndpn*i+5][ndpn*j  ] += Kwu[2][0]; ke[ndpn*i+5][ndpn*j+1] += Kwu[2][1]; ke[ndpn*i+5][ndpn*j+2] += Kwu[2][2];
                
                ke[ndpn*i+3][ndpn*j+3] += Kww[0][0]; ke[ndpn*i+3][ndpn*j+4] += Kww[0][1]; ke[ndpn*i+3][ndpn*j+5] += Kww[0][2];
                ke[ndpn*i+4][ndpn*j+3] += Kww[1][0]; ke[ndpn*i+4][ndpn*j+4] += Kww[1][1]; ke[ndpn*i+4][ndpn*j+5] += Kww[1][2];
                ke[ndpn*i+5][ndpn*j+3] += Kww[2][0]; ke[ndpn*i+5][ndpn*j+4] += Kww[2][1]; ke[ndpn*i+5][ndpn*j+5] += Kww[2][2];
                
                // calculate the kpu matrix
                gp = vec3d(0,0,0);
                for (isol=0; isol<nsol; ++isol) gp += (D[isol]*gradc[isol])*(kappa[isol]/D0[isol]);
                gp = gradp+gp*(R*T);
                wu = vdotTdotv(-gp, dKedE, gradMu[j]);
                ww = vdotTdotv(-gp, dKedE, gradMw[j]);
                for (isol=0; isol<nsol; ++isol) {
                    wu += (((Ke*(D[isol]*gradc[isol])) & gradMu[j])*(J*dkdJ[isol] - kappa[isol])
                           +Ke*(2*kappa[isol]*(gradMu[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
                    + (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradMu[j]))*(-kappa[isol]*R*T/D0[isol]);
                    ww += (((Ke*(D[isol]*gradc[isol])) & gradMw[j])*(J*dkdJ[isol] - kappa[isol])
                           +Ke*(2*kappa[isol]*(gradMw[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
                    + (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradMw[j]))*(-kappa[isol]*R*T/D0[isol]);
                }
                qpu = Phie*gradMu[j];
                qpw = Phie*gradMw[j];
                vec3d kpu = (wu.transpose()*gradMu[i] + qpu*Mu[i])*(detJ*dt);
                vec3d kpw = (ww.transpose()*gradMu[i] + qpw*Mu[i])*(detJ*dt);
                vec3d kqu = (wu.transpose()*gradMw[i] + qpu*Mw[i])*(detJ*dt);
                vec3d kqw = (ww.transpose()*gradMw[i] + qpw*Mw[i])*(detJ*dt);
                ke[ndpn*i+6][ndpn*j  ] += kpu.x; ke[ndpn*i+6][ndpn*j+1] += kpu.y; ke[ndpn*i+6][ndpn*j+2] += kpu.z;
                ke[ndpn*i+6][ndpn*j+3] += kpw.x; ke[ndpn*i+6][ndpn*j+4] += kpw.y; ke[ndpn*i+6][ndpn*j+5] += kpw.z;
                ke[ndpn*i+7][ndpn*j  ] += kqu.x; ke[ndpn*i+7][ndpn*j+1] += kqu.y; ke[ndpn*i+7][ndpn*j+2] += kqu.z;
                ke[ndpn*i+7][ndpn*j+3] += kqw.x; ke[ndpn*i+7][ndpn*j+4] += kqw.y; ke[ndpn*i+7][ndpn*j+5] += kqw.z;
                
                // calculate the kup matrix
                vec3d kup = gradMu[i]*(-Mu[j]*detJ);
                vec3d kuq = gradMu[i]*(-Mw[j]*detJ);
                vec3d kwp = gradMw[i]*(-Mu[j]*detJ);
                vec3d kwq = gradMw[i]*(-Mw[j]*detJ);
                
                ke[ndpn*i  ][ndpn*j+6] += kup.x; ke[ndpn*i  ][ndpn*j+7] += kuq.x;
                ke[ndpn*i+1][ndpn*j+6] += kup.y; ke[ndpn*i+1][ndpn*j+7] += kuq.y;
                ke[ndpn*i+2][ndpn*j+6] += kup.z; ke[ndpn*i+2][ndpn*j+7] += kuq.z;
                
                ke[ndpn*i+3][ndpn*j+6] += kwp.x; ke[ndpn*i+3][ndpn*j+7] += kwq.x;
                ke[ndpn*i+4][ndpn*j+6] += kwp.y; ke[ndpn*i+4][ndpn*j+7] += kwq.y;
                ke[ndpn*i+5][ndpn*j+6] += kwp.z; ke[ndpn*i+5][ndpn*j+7] += kwq.z;
                
                // calculate the kpp matrix
                ke[ndpn*i+6][ndpn*j+6] += (Mu[i]*Mu[j]*Phip - gradMu[i]*(Ke*gradMu[j]))*(detJ*dt);
                ke[ndpn*i+6][ndpn*j+7] += (Mu[i]*Mw[j]*Phip - gradMu[i]*(Ke*gradMw[j]))*(detJ*dt);
                ke[ndpn*i+7][ndpn*j+6] += (Mw[i]*Mu[j]*Phip - gradMw[i]*(Ke*gradMu[j]))*(detJ*dt);
                ke[ndpn*i+7][ndpn*j+7] += (Mw[i]*Mw[j]*Phip - gradMw[i]*(Ke*gradMw[j]))*(detJ*dt);
                
                // calculate kcu matrix data
                jue.zero(); jwe.zero();
                De.zero();
                for (isol=0; isol<nsol; ++isol) {
                    gc[isol] = -gradc[isol]*phiw + w*c[isol]/D0[isol];
                    ju[isol] = ((D[isol]*gc[isol]) & gradMu[j])*(J*dkdJ[isol])
                    + vdotTdotv(gc[isol], dDdE[isol], gradMu[j])*kappa[isol]
                    + (((D[isol]*gradc[isol]) & gradMu[j])*(-phis)
                       +(D[isol]*((gradMu[j]*w)*2) - ((D[isol]*w) & gradMu[j]))*c[isol]/D0[isol]
                       )*kappa[isol]
                    +D[isol]*wu*(kappa[isol]*c[isol]/D0[isol]);
                    jw[isol] = ((D[isol]*gc[isol]) & gradMw[j])*(J*dkdJ[isol])
                    + vdotTdotv(gc[isol], dDdE[isol], gradMw[j])*kappa[isol]
                    + (((D[isol]*gradc[isol]) & gradMw[j])*(-phis)
                       +(D[isol]*((gradMw[j]*w)*2) - ((D[isol]*w) & gradMw[j]))*c[isol]/D0[isol]
                       )*kappa[isol]
                    +D[isol]*ww*(kappa[isol]*c[isol]/D0[isol]);
                    jue += ju[isol]*z[isol];
                    jwe += jw[isol]*z[isol];
                    De += D[isol]*(z[isol]*kappa[isol]*c[isol]/D0[isol]);
                }
                
                for (isol=0; isol<nsol; ++isol) {
                    
                    // calculate the kcu matrix
                    vec3d kcu = ((ju[isol]+jue*penalty).transpose()*gradMu[i])*(detJ*dt);
                    vec3d kcw = ((jw[isol]+jwe*penalty).transpose()*gradMu[i])*(detJ*dt);
                    vec3d kdu = ((ju[isol]+jue*penalty).transpose()*gradMw[i])*(detJ*dt);
                    vec3d kdw = ((jw[isol]+jwe*penalty).transpose()*gradMw[i])*(detJ*dt);
                    ke[ndpn*i+8+2*isol][ndpn*j  ] += kcu.x; ke[ndpn*i+8+2*isol][ndpn*j+1] += kcu.y; ke[ndpn*i+8+2*isol][ndpn*j+2] += kcu.z;
                    ke[ndpn*i+8+2*isol][ndpn*j+3] += kcw.x; ke[ndpn*i+8+2*isol][ndpn*j+4] += kcw.y; ke[ndpn*i+8+2*isol][ndpn*j+5] += kcw.z;
                    ke[ndpn*i+9+2*isol][ndpn*j  ] += kdu.x; ke[ndpn*i+9+2*isol][ndpn*j+1] += kdu.y; ke[ndpn*i+9+2*isol][ndpn*j+2] += kdu.z;
                    ke[ndpn*i+9+2*isol][ndpn*j+3] += kdw.x; ke[ndpn*i+9+2*isol][ndpn*j+4] += kdw.y; ke[ndpn*i+9+2*isol][ndpn*j+5] += kdw.z;
                    
                    // calculate the kcp matrix
                    ke[ndpn*i+8+2*isol][ndpn*j+6] -= (gradMu[i]*(
                                                                 (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                                  +De*penalty)
                                                                 *(Ke*gradMu[j])
                                                                 ))*(detJ*dt);
                    ke[ndpn*i+8+2*isol][ndpn*j+7] -= (gradMu[i]*(
                                                                 (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                                  +De*penalty)
                                                                 *(Ke*gradMw[j])
                                                                 ))*(detJ*dt);
                    ke[ndpn*i+9+2*isol][ndpn*j+6] -= (gradMw[i]*(
                                                                 (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                                  +De*penalty)
                                                                 *(Ke*gradMu[j])
                                                                 ))*(detJ*dt);
                    ke[ndpn*i+9+2*isol][ndpn*j+7] -= (gradMw[i]*(
                                                                 (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                                  +De*penalty)
                                                                 *(Ke*gradMw[j])
                                                                 ))*(detJ*dt);
                    
                    // calculate the kuc matrix
                    sum = 0;
                    for (jsol=0; jsol<nsol; ++jsol)
                        sum += c[jsol]*(dodc[isol]*kappa[jsol]+osmc*dkdc[jsol][isol]);
                    vec3d kuc = (dTdc[isol]*gradMu[i] - gradMu[i]*(R*T*(osmc*kappa[isol]+sum)))*Mu[j]*detJ;
                    vec3d kud = (dTdc[isol]*gradMu[i] - gradMu[i]*(R*T*(osmc*kappa[isol]+sum)))*Mw[j]*detJ;
                    vec3d kwc = (dTdc[isol]*gradMw[i] - gradMw[i]*(R*T*(osmc*kappa[isol]+sum)))*Mu[j]*detJ;
                    vec3d kwd = (dTdc[isol]*gradMw[i] - gradMw[i]*(R*T*(osmc*kappa[isol]+sum)))*Mw[j]*detJ;
                    
                    ke[ndpn*i  ][ndpn*j+8+2*isol] += kuc.x; ke[ndpn*i  ][ndpn*j+9+2*isol] += kud.x;
                    ke[ndpn*i+1][ndpn*j+8+2*isol] += kuc.y; ke[ndpn*i+1][ndpn*j+9+2*isol] += kud.y;
                    ke[ndpn*i+2][ndpn*j+8+2*isol] += kuc.z; ke[ndpn*i+2][ndpn*j+9+2*isol] += kud.z;
                    
                    ke[ndpn*i+3][ndpn*j+8+2*isol] += kwc.x; ke[ndpn*i+3][ndpn*j+9+2*isol] += kwd.x;
                    ke[ndpn*i+4][ndpn*j+8+2*isol] += kwc.y; ke[ndpn*i+4][ndpn*j+9+2*isol] += kwd.y;
                    ke[ndpn*i+5][ndpn*j+8+2*isol] += kwc.z; ke[ndpn*i+5][ndpn*j+9+2*isol] += kwd.z;
                    
                    // calculate the kpc matrix
                    vtmp = vec3d(0,0,0);
                    for (jsol=0; jsol<nsol; ++jsol)
                        vtmp += (D[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol])
                                 +dDdc[jsol][isol]*kappa[jsol])/D0[jsol]*gradc[jsol];
                    wc[isol] = (dKedc[isol]*gp)*(-Mu[j])
                    -Ke*((D[isol]*gradMu[j])*(kappa[isol]/D0[isol])+vtmp*Mu[j])*(R*T);
                    wd[isol] = (dKedc[isol]*gp)*(-Mw[j])
                    -Ke*((D[isol]*gradMw[j])*(kappa[isol]/D0[isol])+vtmp*Mw[j])*(R*T);
                    
                    ke[ndpn*i+6][ndpn*j+8+2*isol] += (gradMu[i]*wc[isol])*(detJ*dt);
                    ke[ndpn*i+6][ndpn*j+9+2*isol] += (gradMu[i]*wd[isol])*(detJ*dt);
                    ke[ndpn*i+7][ndpn*j+8+2*isol] += (gradMw[i]*wc[isol])*(detJ*dt);
                    ke[ndpn*i+7][ndpn*j+9+2*isol] += (gradMw[i]*wd[isol])*(detJ*dt);
                    
                }
                
                // calculate data for the kcc matrix
                jce.assign(nsol, vec3d(0,0,0));
                jde.assign(nsol, vec3d(0,0,0));
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        if (jsol != isol) {
                            jc[isol][jsol] =
                            ((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*Mu[j]
                            +(D[isol]*(w*(-Mu[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            jd[isol][jsol] =
                            ((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*Mw[j]
                            +(D[isol]*(w*(-Mw[j]*dD0dc[isol][jsol]/D0[isol])+wd[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                        }
                        else {
                            jc[isol][jsol] = (D[isol]*(gradMu[j]*(-phiw)+w*(Mu[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*Mu[j]
                            +(D[isol]*(w*(-Mu[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            jd[isol][jsol] = (D[isol]*(gradMw[j]*(-phiw)+w*(Mw[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*Mw[j]
                            +(D[isol]*(w*(-Mw[j]*dD0dc[isol][jsol]/D0[isol])+wd[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                        }
                        jce[jsol] += jc[isol][jsol]*z[isol];
                        jde[jsol] += jd[isol][jsol]*z[isol];
                        
                        // chemical reactions
                        dchatdc[isol][jsol] = 0;
                        for (ireact=0; ireact<nreact; ++ireact)
                            dchatdc[isol][jsol] += m_pMat->GetReaction(ireact)->m_v[isol]
                            *m_pMat->GetReaction(ireact)->Tangent_ReactionSupply_Concentration(mp,jsol);
                        
                        // membrane reactions
                        for (ireact=0; ireact<mreact; ++ireact)
                            dchatdc[isol][jsol] += m_pMat->GetMembraneReaction(ireact)->m_v[isol]
                            *m_pMat->GetMembraneReaction(ireact)->Tangent_ReactionSupply_Concentration(mp,jsol);
                    }
                }
                
                // calculate the kcc matrix
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        ke[ndpn*i+8+2*isol][ndpn*j+8+2*jsol] += (gradMu[i]*(jc[isol][jsol]+jce[jsol]*penalty)
                                                                 + Mu[i]*Mu[j]*phiw*dchatdc[isol][jsol])*(detJ*dt);
                        ke[ndpn*i+8+2*isol][ndpn*j+9+2*jsol] += (gradMu[i]*(jd[isol][jsol]+jde[jsol]*penalty)
                                                                 + Mu[i]*Mw[j]*phiw*dchatdc[isol][jsol])*(detJ*dt);
                        ke[ndpn*i+9+2*isol][ndpn*j+8+2*jsol] += (gradMw[i]*(jc[isol][jsol]+jce[jsol]*penalty)
                                                                 + Mw[i]*Mu[j]*phiw*dchatdc[isol][jsol])*(detJ*dt);
                        ke[ndpn*i+9+2*isol][ndpn*j+9+2*jsol] += (gradMw[i]*(jd[isol][jsol]+jde[jsol]*penalty)
                                                                 + Mw[i]*Mw[j]*phiw*dchatdc[isol][jsol])*(detJ*dt);
                    }
                }
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
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
void FEMultiphasicShellDomain::MembraneReactionFluxes(FEGlobalVector& R)
{
    const int mreact = m_pMat->MembraneReactions();
    if (mreact == 0) return;
    
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // calculate internal force vector
        ElementMembraneReactionFlux(el, fe);
        
        // get the element's LM vector
        UnpackMembraneLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe, true);
    }
}
//-----------------------------------------------------------------------------
//! element external work of flux generated by membrane reactions
void FEMultiphasicShellDomain::ElementMembraneReactionFlux(FEShellElement& el, vector<double>& fe)
{
    double dt = GetFEModel()->GetTime().timeIncrement;

    // get the first material point
    FEMaterialPoint& mp = *el.GetMaterialPoint(0);
    FESolutesMaterialPoint& ps = *mp.ExtractData<FESolutesMaterialPoint>();
    int nse = (int)ps.m_ce.size();
    int nsi = (int)ps.m_ci.size();
    int ndpn = 8 + nse + nsi;

    // nr integration points
    int nint = el.GaussPoints();
    
    // nr of element nodes
    int neln = el.Nodes();
    
    fe.resize(neln*ndpn,0);
    
    // get the element's nodal positions
    vec3d re[FEElement::MAX_NODES], ri[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) {
        FENode& nd = GetFEModel()->GetMesh().Node(el.m_node[j]);
        // nodal positions at front and back surfaces
        re[j] = nd.m_rt;
        ri[j] = nd.st();
    }
    
    double *Mr, *Ms;
    double *M;
    double *w  = el.GaussWeights();
    
    vec3d dxer, dxes, dxet;
    vec3d dxir, dxis, dxit;

    // repeat over integration points
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FESolutesMaterialPoint& ps = *mp.ExtractData<FESolutesMaterialPoint>();
        
        M  = el.H(n);
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        
        // evaluate normal solute flux at integration point
        vector<double> je(nse,0), ji(nsi,0);
        for (int ireact=0; ireact<m_pMat->MembraneReactions(); ++ireact) {
            FEMembraneReaction* react = m_pMat->GetMembraneReaction(ireact);
            FESoluteInterface* psm = react->m_psm;
            double zbar = react->ReactionSupply(mp);
            double zve = 0, zvi = 0;
            for (int k=0; k<nse; ++k) {
                int id = ps.m_ide[k];
                zve += react->m_z[id]*react->m_ve[id];
            }
            for (int k=0; k<nsi; ++k) {
                int id = ps.m_idi[k];
                zvi += react->m_z[id]*react->m_vi[id];
            }
            // Divide fluxes by two because we are integrating over the shell volume
            // but we should only integrate over the shell surface.
            for (int k=0; k<nse; ++k)
                je[k] -= zbar*(react->m_ve[ps.m_ide[k]] + zve)/2;
            for (int k=0; k<nsi; ++k)
                ji[k] -= zbar*(react->m_vi[ps.m_idi[k]] + zvi)/2;
        }
        
        dxer = dxes = vec3d(0,0,0);
        dxir = dxis = vec3d(0,0,0);
        for (int i=0; i<neln; ++i)
        {
            dxer += re[i]*Mr[i];
            dxes += re[i]*Ms[i];
            dxir += ri[i]*Mr[i];
            dxis += ri[i]*Ms[i];
        }
        double ae = (dxer ^ dxes).norm()*w[n]*dt;
        double ai = (dxir ^ dxis).norm()*w[n]*dt;

        for (int i=0; i<neln; ++i)
        {
            for (int k=0; k<nse; ++k)
                fe[ndpn*i+8+k] += M[i]*je[k]*ae;
            for (int k=0; k<nsi; ++k)
                fe[ndpn*i+8+nse+k] += M[i]*ji[k]*ai;
        }
    }

}


//-----------------------------------------------------------------------------
//! calculates the membrane reaction stiffness matrix for this domain
void FEMultiphasicShellDomain::MembraneReactionStiffnessMatrix(FELinearSystem& LS)
{
    const int mreact = m_pMat->MembraneReactions();
    if (mreact == 0) return;
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FEShellElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);

		vector<int> lm;
        UnpackMembraneLM(el, lm);
		ke.SetIndices(lm);
        
        // calculate the element stiffness matrix
        ElementMembraneFluxStiffness(el, ke);
        
        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
//! calculates the element membrane flux stiffness matrix
bool FEMultiphasicShellDomain::ElementMembraneFluxStiffness(FEShellElement& el, matrix& ke)
{
    double dt = GetFEModel()->GetTime().timeIncrement;

    // get the first material point
    FEMaterialPoint& mp = *el.GetMaterialPoint(0);
    FESolutesMaterialPoint& ps = *mp.ExtractData<FESolutesMaterialPoint>();
    int nse = (int)ps.m_ce.size();
    int nsi = (int)ps.m_ci.size();
    int ndpn = 8 + nse + nsi;

    // nr integration points
    int nint = el.GaussPoints();
    
    // nr of element nodes
    int neln = el.Nodes();
    int ndof = ndpn*neln;
    
    // allocate stiffness matrix
    ke.resize(ndof, ndof);
    
    // get the element's nodal coordinates
    vec3d re[FEElement::MAX_NODES], ri[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) {
        FENode& nd = GetFEModel()->GetMesh().Node(el.m_node[j]);
        re[j] = nd.m_rt;
        ri[j] = nd.st();
    }
    
    double *Mr, *Ms;
    double *M;
    double *w  = el.GaussWeights();
    
    vec3d dxir, dxis, dxit;
    vec3d dxer, dxes, dxet;
    
    // repeat over integration points
    ke.zero();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FESolutesMaterialPoint& ps = *mp.ExtractData<FESolutesMaterialPoint>();

        M  = el.H(n);
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        
        // evaluate normal solute flux and its derivatives at integration point
        // take summation over all membrane reactions
        vector<double> je(nse,0), ji(nsi,0);
        vector<double> djedJ(nse,0), djidJ(nsi,0);
        vector<double> djedp(nse,0), djidp(nsi,0);
        vector< vector<double> > djedc(nse,vector<double>(nse,0));
        vector< vector<double> > djidc(nsi,vector<double>(nsi,0));
        for (int ireact=0; ireact<m_pMat->MembraneReactions(); ++ireact) {
            FEMembraneReaction* react = m_pMat->GetMembraneReaction(ireact);
            FESoluteInterface* psm = react->m_psm;
            double zbar = react->ReactionSupply(mp);
            double dzdJ = react->Tangent_ReactionSupply_Strain(mp);
            double dzdpe = react->Tangent_ReactionSupply_Pe(mp);
            double dzdpi = react->Tangent_ReactionSupply_Pi(mp);
            vector<double> dzdce(nse,0), dzdci(nsi,0);
            double zve = 0, zvi = 0;
            for (int k=0; k<nse; ++k) {
                dzdce[k] = react->Tangent_ReactionSupply_Ce(mp, k);
                int id = ps.m_ide[k];
                zve += react->m_z[id]*react->m_ve[id];
            }
            for (int k=0; k<nsi; ++k) {
                dzdci[k] = react->Tangent_ReactionSupply_Ci(mp, k);
                int id = ps.m_idi[k];
                zvi += react->m_z[id]*react->m_vi[id];
            }
            // Divide fluxes by two because we are integrating over the shell volume
            // but we should only integrate over the shell surface.
            for (int k=0; k<nse; ++k) {
                je[k] -= zbar*(react->m_ve[ps.m_ide[k]] + zve)/2;
                djedJ[k] -= dzdJ*(react->m_ve[ps.m_ide[k]] + zve)/2;
                djedp[k] -= dzdpe*(react->m_ve[ps.m_ide[k]] + zve)/2;
                for (int l=0; l<nse; ++l)
                    djedc[k][l] -= dzdce[l]*(react->m_ve[ps.m_ide[k]] + zve)/2;
            }
            for (int k=0; k<nsi; ++k) {
                ji[k] -= zbar*(react->m_vi[ps.m_idi[k]] + zvi)/2;
                djidJ[k] -= dzdJ*(react->m_vi[ps.m_idi[k]] + zvi)/2;
                djidp[k] -= dzdpi*(react->m_vi[ps.m_idi[k]] + zvi)/2;
                for (int l=0; l<nsi; ++l)
                    djidc[k][l] -= dzdci[l]*(react->m_vi[ps.m_idi[k]] + zvi)/2;
            }
        }

        dxer = dxes = vec3d(0,0,0);
        dxir = dxis = vec3d(0,0,0);
        for (int i=0; i<neln; ++i)
        {
            dxer += re[i]*Mr[i];
            dxes += re[i]*Ms[i];
            dxir += ri[i]*Mr[i];
            dxis += ri[i]*Ms[i];
        }
        dxet = dxer ^ dxes;
        dxit = dxir ^ dxis;
        double Je = (dxer ^ dxes).norm();
        double Ji = (dxir ^ dxis).norm();
        vec3d ne = dxet/Je;
        vec3d ni = dxit/Ji;
        
        for (int i=0; i<neln; ++i) {
            int ir = ndpn*i+8;
            for (int j=0; j<neln; ++j)
            {
                int ic = ndpn*j;
                // front face
                vec3d ge = dxes*Mr[j] - dxer*Ms[j];
                mat3d Ae; Ae.skew(ge);
                vec3d ve = Ae*ne;
                for (int k=0; k<nse; ++k) {
                    vec3d kcu = ve*(je[k] + djedJ[k]*Je)*M[i]*w[n]*dt;
                    ke[ir+k][ic  ] -= kcu.x;
                    ke[ir+k][ic+1] -= kcu.y;
                    ke[ir+k][ic+2] -= kcu.z;
                    
                    double kcp = djedp[k]*M[i]*M[j]*Je*w[n]*dt;
                    ke[ir+k][ic+6] -= kcp;
                    
                    for (int l=0; l<nse; ++l) {
                        double kcc = djedc[k][l]*M[i]*M[j]*Je*w[n]*dt;
                        ke[ir+k][ic+8+l] -= kcc;
                    }
                }
                // back face
                vec3d gi = dxis*Mr[j] - dxir*Ms[j];
                mat3d Ai; Ai.skew(gi);
                vec3d vi = Ai*ni;
                for (int k=0; k<nsi; ++k) {
                    vec3d kdw = vi*(ji[k] + djidJ[k]*Ji)*M[i]*w[n]*dt;
                    ke[ir+nse+k][ic+3] -= kdw.x;
                    ke[ir+nse+k][ic+4] -= kdw.y;
                    ke[ir+nse+k][ic+5] -= kdw.z;
                    
                    double kdq = djidp[k]*M[i]*M[j]*Ji*w[n]*dt;
                    ke[ir+nse+k][ic+7] -= kdq;
                    
                    for (int l=0; l<nsi; ++l) {
                        double kdd = djidc[k][l]*M[i]*M[j]*Ji*w[n]*dt;
                        ke[ir+nse+k][ic+8+nse+l] -= kdd;
                    }
                }
            }
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicShellDomain::Update(const FETimeInfo& tp)
{
	FESSIShellDomain::Update(tp);

    bool berr = false;
    int NE = (int) m_Elem.size();
#pragma omp parallel for shared(NE, berr)
    for (int i=0; i<NE; ++i)
    {
        try
        {
            UpdateElementStress(i, tp);
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
void FEMultiphasicShellDomain::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double dt = tp.timeIncrement;
    
    int j, k, n;
    int nint, neln;
    double* gw;
    vec3d r0[FEElement::MAX_NODES];
    vec3d rt[FEElement::MAX_NODES];
    double pn[FEElement::MAX_NODES], qn[FEElement::MAX_NODES];
    DOFS& fedofs = GetFEModel()->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    int MAX_DDOFS = fedofs.GetVariableSize("shell concentration");

    FEMesh& mesh = *m_pMesh;
    
    // get the multiphasic material
    FEMultiphasic* pmb = m_pMat;
    const int nsol = (int)pmb->Solutes();
    vector<int> sid(nsol);
    for (j=0; j<nsol; ++j) sid[j] = pmb->GetSolute(j)->GetSoluteDOF();
    
    // get the shell element
    FEShellElement& el = m_Elem[iel];
    
    // get the number of integration points
    nint = el.GaussPoints();
    
    // get the number of nodes
    neln = el.Nodes();
    vector< vector<double> > cn(MAX_CDOFS, vector<double>(neln));
    vector< vector<double> > dn(MAX_DDOFS, vector<double>(neln));

    // get the integration weights
    gw = el.GaussWeights();
    
    const double *M, *Mr, *Ms;
    
    // get the nodal data
    for (j=0; j<neln; ++j)
    {
        r0[j] = mesh.Node(el.m_node[j]).m_r0;
        rt[j] = mesh.Node(el.m_node[j]).m_rt;
        pn[j] = mesh.Node(el.m_node[j]).get(m_dofP);
        qn[j] = mesh.Node(el.m_node[j]).get(m_dofQ);
        for (k=0; k<MAX_CDOFS; ++k)
            cn[k][j] = mesh.Node(el.m_node[j]).get(m_dofC + k);
        for (k=0; k<MAX_DDOFS; ++k)
            dn[k][j] = mesh.Node(el.m_node[j]).get(m_dofD + k);
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
        
        // update membrane reaction data if needed
        if (m_pMat->MembraneReactions()) {
            int nse = (int)spt.m_ce.size();
            int nsi = (int)spt.m_ci.size();
            M = el.H(n);
            Mr = el.Hr(n);
            Ms = el.Hs(n);
            double pe = 0, pi = 0;
            vector<double> ce(nse,0), ci(nsi,0);
            vec3d dxr(0,0,0), dxs(0,0,0);
            vec3d dXr(0,0,0), dXs(0,0,0);
            for (int j=0; j<neln; ++j) {
                dxr += rt[j]*Mr[j];
                dxs += rt[j]*Ms[j];
                dXr += r0[j]*Mr[j];
                dXs += r0[j]*Ms[j];
                pe += pn[j]*M[j];
                pi += qn[j]*M[j];
                for (int k=0; k<nse; ++k)
                    ce[k] += cn[spt.m_ide[k]][j]*M[j];
                for (int k=0; k<nsi; ++k)
                    ci[k] += dn[spt.m_idi[k]][j]*M[j];
            }
            spt.m_strain = (dxr ^ dxs).norm()/(dXr ^ dXs).norm();
            spt.m_pe = pe;
            spt.m_pi = pi;
            spt.m_ce = ce;
            spt.m_ci = ci;
        }
        
        for (k=0; k<nsol; ++k) {
            // evaluate effective solute concentrations at gauss-point
            spt.m_c[k] = evaluate(el, cn[sid[k]], dn[sid[k]], n);
            // calculate the gradient of c at gauss-point
            spt.m_gradc[k] = gradient(el, cn[sid[k]], dn[sid[k]], n);
        }
        
        // update SBM referential densities
        pmb->UpdateSolidBoundMolecules(mp);
        
        // evaluate referential solid volume fraction
        ppt.m_phi0t = pmb->SolidReferentialVolumeFraction(mp);
        if (m_breset) ppt.m_phi0 = ppt.m_phi0t;
        
        // evaluate fluid pressure at gauss-point
        ppt.m_p = evaluate(el, pn, qn, n);
        
        // calculate the gradient of p at gauss-point
        ppt.m_gradp = gradient(el, pn, qn, n);
        
        // update the fluid and solute fluxes
        // and evaluate the actual fluid pressure and solute concentration
        ppt.m_w = pmb->FluidFlux(mp);
        spt.m_psi = pmb->ElectricPotential(mp);
        for (k=0; k<nsol; ++k) {
            spt.m_ca[k] = pmb->Concentration(mp,k);
            spt.m_j[k] = pmb->SoluteFlux(mp,k);
        }
        ppt.m_pa = pmb->Pressure(mp);
        spt.m_cF = pmb->FixedChargeDensity(mp);
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
        
        // update membrane reaction element data
        for (int j=0; j<m_pMat->MembraneReactions(); ++j)
            pmb->GetMembraneReaction(j)->UpdateElementData(mp);
        
    }
    if (m_breset) m_breset = false;
}

//-----------------------------------------------------------------------------
// Extract the solute DOFs for the solid elements on either side of the shell domain
void FEMultiphasicShellDomain::UpdateShellMPData(int iel)
{
    if (m_pMat->MembraneReactions() == 0) return;
    
    static bool bfirst = true;
    
    if (bfirst) {
        FEMesh& mesh = *GetMesh();

        // get the shell element
        FEShellElement& el = m_Elem[iel];
        FEElement* si = mesh.FindElementFromID(el.m_elem[0]);
        FEElement* se = mesh.FindElementFromID(el.m_elem[1]);
        if ((si == nullptr) || (se == nullptr)) return;
        FEMultiphasic* mpi = dynamic_cast<FEMultiphasic*>(GetFEModel()->GetMaterial(si->GetMatID()));
        FEMultiphasic* mpe = dynamic_cast<FEMultiphasic*>(GetFEModel()->GetMaterial(se->GetMatID()));
        if ((mpi == nullptr) || (mpe == nullptr)) return;

        vector<int> idi(mpi->Solutes(),-1);
        vector<int> ide(mpe->Solutes(),-1);
        for (int i=0; i<mpi->Solutes(); ++i) idi[i] = mpi->GetSolute(i)->GetSoluteDOF();
        for (int i=0; i<mpe->Solutes(); ++i) ide[i] = mpe->GetSolute(i)->GetSoluteDOF();

        // get the number of integration points
        int nint = el.GaussPoints();
        
        // loop over the integration points
        for (int n = 0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            
            ps.m_ide = ide;
            ps.m_idi = idi;
            ps.m_ce.resize(ide.size());
            ps.m_ci.resize(idi.size());
        }
        
        bfirst = false;
    }
}

