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

#include "stdafx.h"
#include "FEMultiphasicFSIDomain3D.h"
#include <FECore/log.h>
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioMultiphasicFSI.h"
#include "FEFluidFSI.h"
#include "FEBiphasicFSI.h"
#include <FECore/FELinearSystem.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEMultiphasicFSIDomain3D::FEMultiphasicFSIDomain3D(FEModel* pfem) : FESolidDomain(pfem), FEMultiphasicFSIDomain(pfem), m_dofU(pfem), m_dofV(pfem), m_dofW(pfem), m_dofAW(pfem), m_dofSU(pfem), m_dofR(pfem), m_dof(pfem)
{
    m_pMat = 0;
    m_btrans = true;
    m_sseps = 0;
    
    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofU.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::DISPLACEMENT));
        m_dofV.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::VELOCITY));
        m_dofW.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::RELATIVE_FLUID_VELOCITY));
        m_dofAW.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::RELATIVE_FLUID_ACCELERATION));
        m_dofSU.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::SHELL_DISPLACEMENT));
        m_dofR.AddVariable(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::RIGID_ROTATION));
        m_dofEF = pfem->GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_DILATATION), 0);
        m_dofAEF = pfem->GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_DILATATION_TDERIV), 0);
        m_dofC = pfem->GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_CONCENTRATION), 0);
        m_dofAC = pfem->GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_CONCENTRATION_TDERIV), 0);
    }
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEMultiphasicFSIDomain3D& FEMultiphasicFSIDomain3D::operator = (FEMultiphasicFSIDomain3D& d)
{
    m_Elem = d.m_Elem;
    m_pMesh = d.m_pMesh;
    return (*this);
}

//-----------------------------------------------------------------------------
// get the total dof
const FEDofList& FEMultiphasicFSIDomain3D::GetDOFList() const
{
    return m_dof;
}

//-----------------------------------------------------------------------------
//! Assign material
void FEMultiphasicFSIDomain3D::SetMaterial(FEMaterial* pmat)
{
    FEDomain::SetMaterial(pmat);
    if (pmat)
    {
        m_pMat = dynamic_cast<FEMultiphasicFSI*>(pmat);
        assert(m_pMat);
    }
    else m_pMat = 0;
}

//-----------------------------------------------------------------------------
bool FEMultiphasicFSIDomain3D::Init()
{
    // initialize base class
    if (FESolidDomain::Init() == false) return false;
    
    //TODO Initialize body force like biphasic solver? maybe not necessary.
    const int nsol = m_pMat->Solutes();
    
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
            FEBiphasicFSIMaterialPoint& pb = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
            FEMultiphasicFSIMaterialPoint& ps = *(mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
            
            pb.m_phi0 = m_pMat->SolidReferentialVolumeFraction(mp);
            ps.m_cF = m_pMat->FixedChargeDensity(mp);
        }
    }
    
    // set the active degrees of freedom list
    FEDofList dofs(GetFEModel());
    for (int i=0; i<nsol; ++i)
    {
        int m = m_pMat->GetSolute(i)->GetSoluteDOF();
        dofs.AddDof(m_dofC + m);
    }
    m_dof = dofs;
    return true;
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::Activate()
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
            }
            node.set_active(m_dofW[0]);
            node.set_active(m_dofW[1]);
            node.set_active(m_dofW[2]);
            node.set_active(m_dofEF);
            for (int isol=0; isol<nsol; ++isol)
                node.set_active(m_dofC + m_pMat->GetSolute(isol)->GetSoluteDOF());
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::InitMaterialPoints()
{
    const int nsol = m_pMat->Solutes();
    FEMesh& m = *GetMesh();
    
    const int NE = FEElement::MAX_NODES;
    double ef[NE];
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
        for (int i = 0; i<neln; ++i)
        {
            FENode& ni = m.Node(el.m_node[i]);
            ef[i] = ni.get(m_dofEF);
            for (int isol = 0; isol<nsol; ++isol)
                c0[isol][i] = ni.get(m_dofC + sid[isol]);
        }
        
        // get the number of integration points
        int nint = el.GaussPoints();
        
        // loop over the integration points
        for (int n = 0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEFluidMaterialPoint& ft = *(mp.ExtractData<FEFluidMaterialPoint>());
            FEFSIMaterialPoint& fs = *(mp.ExtractData<FEFSIMaterialPoint>());
            FEBiphasicFSIMaterialPoint& pt = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
            FEMultiphasicFSIMaterialPoint& ps = *(mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
            
            // initialize effective fluid pressure, its gradient, and fluid flux
            ft.m_ef = el.Evaluate(ef, n);
            ft.m_gradef = gradient(el, ef, n);
            ps.m_pe = m_pMat->Fluid()->Pressure(mp);
            
            // initialize solutes
            ps.m_nsol = nsol;
            
            // initialize effective solute concentrations
            for (int isol = 0; isol<nsol; ++isol) {
                ps.m_c[isol] = el.Evaluate(c0[isol], n);
                ps.m_gradc[isol] = gradient(el, c0[isol], n);
            }
            
            ps.m_psi = m_pMat->ElectricPotential(mp);
            for (int isol = 0; isol<nsol; ++isol) {
                ps.m_ca[isol] = m_pMat->ConcentrationActual(mp, isol);
                ps.m_j[isol] = m_pMat->SoluteFlux(mp, isol);
            }
            ft.m_pf = m_pMat->PressureActual(mp);
            
            // initialize referential solid volume fraction
            pt.m_phi0 = m_pMat->SolidReferentialVolumeFraction(mp);
            
            // calculate FCD, current and stress
            ps.m_cF = m_pMat->FixedChargeDensity(mp);
            ps.m_Ie = m_pMat->CurrentDensity(mp);
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::Reset()
{
    // reset base class data
    FESolidDomain::Reset();
    
    const int nsol = m_pMat->Solutes();
    
    // initialize all element data
    ForEachMaterialPoint([=](FEMaterialPoint& mp) {
        FEBiphasicFSIMaterialPoint& pt = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
        FEMultiphasicFSIMaterialPoint& ps = *(mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
        
        // initialize referential solid volume fraction
        pt.m_phi0 = m_pMat->m_phi0(mp);
        
        // initialize solutes
        ps.m_nsol = nsol;
        ps.m_c.assign(nsol,0);
        ps.m_ca.assign(nsol,0);
        ps.m_cdot.assign(nsol,0);
        ps.m_gradc.assign(nsol,vec3d(0,0,0));
        ps.m_j.assign(nsol,vec3d(0,0,0));
        ps.m_k.assign(nsol, 0);
        ps.m_dkdJ.assign(nsol, 0);
        ps.m_dkdc.resize(nsol, vector<double>(nsol,0));
    });
    
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
            for (int j=0; j<m_pMat->Reactions(); ++j)
                m_pMat->GetReaction(j)->ResetElementData(mp);
        }
    }
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEMultiphasicFSIDomain3D::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    const int NE = FEElement::MAX_NODES;
    vec3d x0[NE], xt[NE], r0, rt, v;
    FEMesh& m = *GetMesh();
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FESolidElement& el = m_Elem[i];
        if (el.isActive())
        {
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
                FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
                FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
                
                mp.m_r0 = r0;
                mp.m_rt = rt;
                et.m_Wp = et.m_Wt;
                
                if ((pt.m_ef <= -1) || (et.m_J <= 0)) {
                    throw NegativeJacobianDetected();
                }
                
                // reset chemical reaction element data
                for (int j=0; j<m_pMat->Reactions(); ++j)
                    m_pMat->GetReaction(j)->InitializeElementData(mp);
                
                mp.Update(timeInfo);
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data.
void FEMultiphasicFSIDomain3D::UnpackLM(FEElement& el, vector<int>& lm)
{
    int N = el.Nodes();
    int nsol = m_pMat->Solutes();
    int ndpn = 7+nsol;
    lm.resize(N*(3+ndpn));
    for (int i=0; i<N; ++i)
    {
        FENode& node = m_pMesh->Node(el.m_node[i]);
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[ndpn*i  ] = id[m_dofU[0]];
        lm[ndpn*i+1] = id[m_dofU[1]];
        lm[ndpn*i+2] = id[m_dofU[2]];
        lm[ndpn*i+3] = id[m_dofW[0]];
        lm[ndpn*i+4] = id[m_dofW[1]];
        lm[ndpn*i+5] = id[m_dofW[2]];
        lm[ndpn*i+6] = id[m_dofEF];
        for (int isol = 0; isol < nsol; ++isol)
        {
            //int m = m_pMat->GetSolute(i)->GetSoluteDOF();
            lm[ndpn*i+7+isol] = id[m_dofC + m_pMat->GetSolute(isol)->GetSoluteDOF()];
        }
        
        // rigid rotational dofs
        lm[ndpn*N + 3*i  ] = id[m_dofR[0]];
        lm[ndpn*N + 3*i+1] = id[m_dofR[1]];
        lm[ndpn*N + 3*i+2] = id[m_dofR[2]];
    }
    
    // substitute interface dofs for solid-shell interfaces
    FESolidElement& sel = static_cast<FESolidElement&>(el);
    for (int i = 0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = m_pMesh->Node(el.m_node[i]);
            vector<int>& id = node.m_ID;
            
            // first the displacement dofs
            lm[ndpn*i  ] = id[m_dofSU[0]];
            lm[ndpn*i+1] = id[m_dofSU[1]];
            lm[ndpn*i+2] = id[m_dofSU[2]];
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::InternalForces(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
    
    int nsol = m_pMat->Solutes();
    int ndpn = 7+nsol;
    
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        if (el.isActive()) {
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
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEMultiphasicFSIDomain3D::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJ;
    
    mat3ds sv, se, pi;
    vec3d gradp, divTe;
    mat3ds km1;
    
    const double *H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    int ndpn = 7+nsol;
    const int nreact = m_pMat->Reactions();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double*    gw = el.GaussWeights();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEElasticMaterialPoint& et = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEFSIMaterialPoint& ft = *(mp.ExtractData<FEFSIMaterialPoint>());
        FEBiphasicFSIMaterialPoint& bt = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
        FEMultiphasicFSIMaterialPoint& mt = *(mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
        double Jf =  1 + pt.m_ef;
        
        // calculate the jacobian
        detJ = invjact(el, Ji, n, tp.alphaf)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        // get the viscous stress tensor for this integration point
        sv = m_pMat->Fluid()->GetViscous()->Stress(mp);
        se = m_pMat->Solid()->Stress(mp);
        double pa = m_pMat->PressureActual(mp);
        // get the gradient of the elastic pressure
        gradp = pt.m_gradef*m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        // get inverse of permeability tensor
        km1 = m_pMat->InvPermeability(mp);
        //get pI
        pi = mat3dd(m_pMat->Fluid()->Pressure(mp));
        
        // Miscellaneous constants
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double penalty = m_pMat->m_penalty;
        double dms = m_pMat->m_diffMtmSupp;
        
        // evaluate the chat
        vector<double> chat(nsol,0);
        double phiwhat = 0;
        
        // chemical reactions
        for (i=0; i<nreact; ++i) {
            FEChemicalReaction* pri = m_pMat->GetReaction(i);
            double zhat = pri->ReactionSupply(mp);
            phiwhat += pri->m_Vbar*zhat;
            for (int isol=0; isol<nsol; ++isol)
            {
                chat[isol] += zhat*pri->m_v[isol];
            }
        }
        
        vector<double> M(nsol);
        vector<int> z(nsol);
        vec3d je(0,0,0);
        vector<double> d0(nsol);
        vector<mat3ds> dm1(nsol);
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            M[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            dm1[isol] = m_pMat->InvDiffusivity(mp, isol);
            je += mt.m_j[isol]*z[isol];
        }
        
        vector<double> dkdt(nsol,0);
        for (int isol=0; isol<nsol; ++isol)
        {
            dkdt[isol] = mt.m_dkdJ[isol]*ft.m_Jdot;
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                dkdt[isol] += mt.m_dkdc[isol][jsol]*mt.m_cdot[jsol];
            }
        }
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
        {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        // Jfdot wrt fluid
        double phif = m_pMat->Porosity(mp);
        double phis = m_pMat->SolidVolumeFrac(mp);
        vec3d gradphif = m_pMat->gradPorosity(mp);
        vec3d gradphifphis = m_pMat->gradPhifPhis(mp);
        double dJfdotf = pt.m_efdot + pt.m_gradef*ft.m_w/phif;
        // Jsdot/Js
        double dJsoJ = ft.m_Jdot/et.m_J;
        double phifdot = dJsoJ*phis;
        mat3dd I = mat3dd(1.0);
        
         //Old Nat BC
        mat3ds Tmix = se - sv*phis;
        for (int isol = 0; isol < nsol; ++isol)
        {
            Tmix += -mat3dd(1.0)*R*T*osmc*mt.m_ca[isol];
        }
        
        
        //New Nat BC
        //mat3ds Tmix = -mat3dd(1.0)*pa + se - sv*phis;
        
        for (i=0; i<neln; ++i)
        {
            vec3d fs = (Tmix*gradN[i] + (sv*gradphifphis*phis*phis/phif - km1*ft.m_w)*H[i])*detJ; //Old Nat BC
            //vec3d fs = (Tmix*gradN[i] + (-gradp + sv*gradphifphis*phis*phis/phif - km1*ft.m_w)*H[i])*detJ; //New Nat BC
            vec3d ff = (sv*gradN[i] + (gradp + km1*ft.m_w - sv*gradphif/phif)*H[i])*detJ;
            double fJ = (H[i]*(dJfdotf*phif/Jf - dJsoJ + phif*phiwhat) + gradN[i]*ft.m_w)*detJ;
            
            for (int isol=0; isol<nsol; ++isol)
            {
                fs += (-mt.m_gradc[isol]*mt.m_k[isol]*phif*R*T - (dm1[isol] - I/phif/d0[isol])*mt.m_j[isol]*R*T - ft.m_w*R*T*phis/phif*mt.m_k[isol]*mt.m_c[isol]/d0[isol])*H[i]*detJ*dms;
                ff += (ft.m_w*R*T/phif*mt.m_k[isol]*mt.m_c[isol]/d0[isol] - mt.m_j[isol]*R*T/phif/d0[isol])*H[i]*dms*detJ;
            }
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[ndpn*i  ] -= fs.x;
            fe[ndpn*i+1] -= fs.y;
            fe[ndpn*i+2] -= fs.z;
            fe[ndpn*i+3] -= ff.x;
            fe[ndpn*i+4] -= ff.y;
            fe[ndpn*i+5] -= ff.z;
            fe[ndpn*i+6] -= fJ;
            
            for (int isol=0; isol<nsol; ++isol)
            {
                double fc = ((mt.m_j[isol]+je*penalty)*gradN[i] + H[i]*(phif*chat[isol] - (ft.m_Jdot*phif*mt.m_k[isol]*mt.m_c[isol] + phifdot*et.m_J*mt.m_k[isol]*mt.m_c[isol] + dkdt[isol]*et.m_J*phif*mt.m_c[isol] + mt.m_cdot[isol]*et.m_J*phif*mt.m_k[isol])/et.m_J))*detJ;
                fe[ndpn*i+7+isol] -= fc;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
    
    int nsol = m_pMat->Solutes();
    int ndpn = 7+nsol;
    for (int i=0; i<NE; ++i)
    {
        // get the element
        FESolidElement& el = m_Elem[i];
        
        if (el.isActive()) {
            vector<double> fe;
            vector<int> lm;
            
            // get the element force vector and initialize it to zero
            int ndof = ndpn*el.Nodes();
            fe.assign(ndof, 0);
            
            // apply body forces
            ElementBodyForce(BF, el, fe);
            
            // get the element's LM vector
            UnpackLM(el, lm);
            
            // assemble element 'fe'-vector into global R vector
            R.Assemble(el.m_node, lm, fe);
        }
    }
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEMultiphasicFSIDomain3D::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    // jacobian
    double Ji[3][3], detJ;
    double *H, *Gr, *Gs, *Gt;
    double* gw = el.GaussWeights();
    double cf;
    vec3d ff, f;
    
    int neln = el.Nodes();
    int nsol = m_pMat->Solutes();
    int ndpn = 7+nsol;
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    // nodal coordinates
    vec3d r0[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
        r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
        FEMultiphasicFSIMaterialPoint& mt = *mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
        
        double densTs = m_pMat->TrueSolidDensity(mp);
        double densTf = m_pMat->TrueFluidDensity(mp);
        
        pt.m_r0 = el.Evaluate(r0, n);
        
        detJ = invjact(el, Ji, n, tp.alphaf)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double dms = m_pMat->m_diffMtmSupp;
        double penalty = m_pMat->m_penalty;
        vector<double> M(nsol);
        vector<int> z(nsol);
        vector<double> d0(nsol);
        vector<mat3ds> D(nsol);
        
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            M[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            D[isol] = m_pMat->Diffusivity(mp, isol);
        }
        
        // evaluate spatial gradient of shape functions
        for (int i=0; i<neln; ++i)
        {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        // get the force
        double phis = m_pMat->SolidVolumeFrac(mp);
        double phif = m_pMat->Porosity(mp);
        mat3dd I = mat3dd(1.0);
        
        for (int i=0; i<neln; ++i)
        {
            f = BF.force(mp)*((-densTf+densTs)*phis*detJ*H[i]);
            ff = BF.force(mp)*(densTf*detJ*H[i]);
            
            for (int isol = 0; isol < nsol; ++isol)
            {
                f += (I*phif-D[isol]/d0[isol])*BF.force(mp)*M[isol]*mt.m_k[isol]*mt.m_c[isol]*H[i]*detJ*dms;
                ff+= D[isol]*BF.force(mp)/d0[isol]*M[isol]*mt.m_k[isol]*mt.m_c[isol]*H[i]*detJ*dms;
            }
            
            fe[ndpn*i+0] -= f.x;
            fe[ndpn*i+1] -= f.y;
            fe[ndpn*i+2] -= f.z;
            fe[ndpn*i+3] -= ff.x;
            fe[ndpn*i+4] -= ff.y;
            fe[ndpn*i+5] -= ff.z;
            
            for (int isol=0; isol<nsol; ++isol)
            {
                double fc = -gradN[i]*(D[isol]*BF.force(mp))*mt.m_k[isol]*M[isol]*mt.m_c[isol]*phif/(R*T)*detJ;
                for(int jsol = 0; jsol<nsol; ++jsol)
                {
                    fc += -gradN[i]*(D[jsol]*BF.force(mp))*z[jsol]*penalty*mt.m_k[jsol]*M[jsol]*mt.m_c[jsol]*phif/(R*T)*detJ;
                }
                fe[ndpn*i+7+isol] -= fc;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
//! For now, we assume that the body force is constant
void FEMultiphasicFSIDomain3D::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int neln = el.Nodes();
    
    // jacobian
    double Ji[3][3], detJ;
    double *H, *Gr, *Gs, *Gt;
    double* gw = el.GaussWeights();
    vec3d f, ff, fs, ffs, kwJ, kuJ;
    mat3d Kwu, Kuu;
    int nsol = m_pMat->Solutes();
    int ndpn = 7+nsol;
    
    const int nreact = m_pMat->Reactions();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
        FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
        FEMultiphasicFSIMaterialPoint& mt = *(mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
        double Jf = 1 + pt.m_ef;
        
        // calculate the jacobian
        detJ = invjact(el, Ji, n, tp.alphaf)*gw[n]*tp.alphaf;
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        
        double densf = m_pMat->FluidDensity(mp);
        double densTf = m_pMat->TrueFluidDensity(mp);
        double denss = m_pMat->SolidDensity(mp);
        double densTs = m_pMat->TrueSolidDensity(mp);
        
        // get the force
        ff = BF.force(mp)*(densTf*detJ);
        f = BF.force(mp)*(densf*detJ);
        fs = BF.force(mp)*(densTs*detJ);
        double phif = m_pMat->Porosity(mp);
        double phis = m_pMat->SolidVolumeFrac(mp);
        ffs = BF.force(mp)*((phis/phif*(densf+denss) - densTf*phis + denss)*detJ);
        
        double dens = m_pMat->Fluid()->Density(mp);
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double dms = m_pMat->m_diffMtmSupp;
        double penalty = m_pMat->m_penalty;
        vector<double> M(nsol);
        vector<int> z(nsol);
        vector<double> d0(nsol);
        vector<vector<double>> d0p(nsol, vector<double>(nsol));
        vector<mat3ds> D(nsol);
        vector<tens4dmm> dDdE(nsol);
        vector<vector<mat3ds>> dDdc(nsol, vector<mat3ds>(nsol));
        
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            M[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            D[isol] = m_pMat->Diffusivity(mp, isol);
            dDdE[isol] = m_pMat->Diffusivity_Tangent_Strain(mp, isol);
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                d0p[isol][jsol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, jsol);
                dDdc[isol][jsol] = m_pMat->Diffusivity_Tangent_Concentration(mp, isol, jsol);
            }
        }
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (int i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        for (int i=0, i7=0; i<neln; ++i, i7 += ndpn) {
            for (int j=0, j7 = 0; j<neln; ++j, j7 += ndpn)
            {
                kwJ = ff*(-H[i]*H[j]/Jf);
                kuJ = ff*(H[i]*H[j]*phis/Jf);
                Kwu = (ff & gradN[j])*H[i];
                Kuu = mat3d(0.0);
                
                for (int isol = 0; isol < nsol; ++isol)
                {
                    Kwu += (((D[isol]*BF.force(mp))&gradN[j])*M[isol]*mt.m_c[isol]*mt.m_dkdJ[isol]/d0[isol]*H[i]*et.m_J + (vdotTdotv(BF.force(mp), dDdE[isol], gradN[j]) + mat3dd(1.0)*((D[isol]*BF.force(mp))*gradN[j]) + D[isol]*(gradN[j] & BF.force(mp)))*M[isol]*mt.m_k[isol]*mt.m_c[isol]/d0[isol]*H[i])*detJ*dms;
                    Kuu += ((((mat3dd(1.0)*phif-D[isol]/d0[isol])*BF.force(mp))&gradN[j])*M[isol]*mt.m_c[isol]*mt.m_dkdJ[isol]*H[i]*et.m_J + ((BF.force(mp)&gradN[j]) - mat3dd(1.0)*((D[isol]*BF.force(mp))*gradN[j])/d0[isol])*M[isol]*mt.m_c[isol]*mt.m_k[isol]*H[i] - (vdotTdotv(BF.force(mp), dDdE[isol], gradN[j]) + D[isol]*(gradN[j] & BF.force(mp)))*M[isol]*mt.m_k[isol]*mt.m_c[isol]/d0[isol]*H[i])*detJ*dms;
                }
                
                ke[i7+0][j7+0] += Kuu(0,0); ke[i7+0][j7+1] += Kuu(0,1); ke[i7+0][j7+2] += Kuu(0,2);
                ke[i7+1][j7+0] += Kuu(1,0); ke[i7+1][j7+1] += Kuu(1,1); ke[i7+1][j7+2] += Kuu(1,2);
                ke[i7+2][j7+0] += Kuu(2,0); ke[i7+2][j7+1] += Kuu(2,1); ke[i7+2][j7+2] += Kuu(2,2);
                ke[i7+3][j7+0] += Kwu(0,0); ke[i7+3][j7+1] += Kwu(0,1); ke[i7+3][j7+2] += Kwu(0,2);
                ke[i7+4][j7+0] += Kwu(1,0); ke[i7+4][j7+1] += Kwu(1,1); ke[i7+4][j7+2] += Kwu(1,2);
                ke[i7+5][j7+0] += Kwu(2,0); ke[i7+5][j7+1] += Kwu(2,1); ke[i7+5][j7+2] += Kwu(2,2);
                ke[i7+0][j7+6] += kuJ.x;
                ke[i7+1][j7+6] += kuJ.y;
                ke[i7+2][j7+6] += kuJ.z;
                ke[i7+3][j7+6] += kwJ.x;
                ke[i7+4][j7+6] += kwJ.y;
                ke[i7+5][j7+6] += kwJ.z;
                
                for (int isol = 0; isol<nsol; ++isol)
                {
                    vec3d kwc = vec3d(0.0);
                    vec3d kuc = vec3d(0.0);
                    vec3d kcu = ((gradN[i]&gradN[j])*D[isol]*BF.force(mp)*M[isol]*mt.m_k[isol]*mt.m_c[isol]/R/T*phif - gradN[j]*(gradN[i]*(D[isol]*BF.force(mp)))*M[isol]*mt.m_k[isol]*mt.m_c[isol]/R/T*phis -  ((gradN[j]&BF.force(mp))*D[isol]*et.m_J*mt.m_dkdJ[isol] + mat3dd(1.0)*(D[isol]*BF.force(mp)*gradN[j])*mt.m_k[isol] + vdotTdotv(BF.force(mp), dDdE[isol], gradN[j]).transpose()*mt.m_k[isol] + (BF.force(mp)&gradN[j])*D[isol]*mt.m_k[isol]) *gradN[i]*M[isol]*mt.m_c[isol]/R/T*phif)*detJ*dms;
                    for(int jsol=0; jsol<nsol; ++jsol)
                    {
                        double kcc = 0;
                        kcu += ((gradN[i]&gradN[j])*D[jsol]*BF.force(mp)*M[jsol]*z[jsol]*penalty*mt.m_k[jsol]*mt.m_c[jsol]/R/T*phif - gradN[j]*(gradN[i]*(D[jsol]*BF.force(mp)))*M[jsol]*z[jsol]*penalty*mt.m_k[jsol]*mt.m_c[jsol]/R/T*phis - ((gradN[j]&BF.force(mp))*D[jsol]*et.m_J*mt.m_dkdJ[jsol] + mat3dd(1.0)*(D[jsol]*BF.force(mp)*gradN[j])*mt.m_k[jsol] + vdotTdotv(BF.force(mp), dDdE[jsol], gradN[j]).transpose()*mt.m_k[jsol] + (BF.force(mp)&gradN[j])*D[jsol]*mt.m_k[jsol]) *gradN[i]*M[jsol]*z[jsol]*penalty*mt.m_c[jsol]/R/T*phif)*detJ*dms;
                        if (isol == jsol)
                        {
                            kwc += (D[isol]*mt.m_dkdc[isol][isol]*mt.m_c[isol]/d0[isol] + D[isol]*mt.m_k[isol]/d0[isol] - D[isol]*mt.m_k[isol]*mt.m_c[isol]*d0p[isol][isol]/d0[isol]/d0[isol] + dDdc[isol][isol]*mt.m_k[isol]*mt.m_c[isol]/d0[isol])*BF.force(mp)*M[isol]*H[i]*H[j]*detJ*dms;
                            kuc += ((mat3dd(1.0)*phif - D[isol]/d0[isol])*BF.force(mp)*(mt.m_dkdc[isol][isol]*mt.m_c[isol] + mt.m_k[isol])*M[isol]*H[i]*H[j] + (D[isol]*d0p[isol][isol]/d0[isol]/d0[isol] - dDdc[isol][isol]/d0[isol])*BF.force(mp)*M[isol]*mt.m_k[isol]*mt.m_c[isol]*H[i]*H[j])*dms*detJ;
                            kcc = -(D[isol]*mt.m_dkdc[isol][isol]*mt.m_c[isol] + D[isol]*mt.m_k[isol] + dDdc[isol][isol]*mt.m_k[isol]*mt.m_c[isol])*BF.force(mp)*gradN[i]*M[isol]/R/T*H[j]*phif*detJ;
                            for (int ksol=0; ksol<nsol; ++ksol)
                            {
                                if(isol==ksol)
                                {
                                    kcc += -(D[isol]*mt.m_dkdc[isol][isol]*mt.m_c[isol] + D[isol]*mt.m_k[isol] + dDdc[isol][isol]*mt.m_k[isol]*mt.m_c[isol])*BF.force(mp)*gradN[i]*M[isol]*z[isol]*penalty/R/T*H[j]*phif*detJ;
                                }
                                else
                                {
                                    kcc += -(D[ksol]*mt.m_dkdc[ksol][isol]*mt.m_c[ksol] + dDdc[ksol][isol]*mt.m_k[ksol]*mt.m_c[ksol])*BF.force(mp)*gradN[i]*M[ksol]*z[ksol]*penalty/R/T*H[j]*phif*detJ;
                                }
                            }
                        }
                        else
                        {
                            kwc += (D[jsol]*mt.m_dkdc[jsol][isol]*mt.m_c[jsol]/d0[jsol] - D[jsol]*mt.m_k[jsol]*mt.m_c[jsol]*d0p[jsol][isol]/d0[jsol]/d0[jsol] + dDdc[jsol][isol]*mt.m_k[jsol]*mt.m_c[jsol]/d0[jsol])*BF.force(mp)*M[jsol]*H[i]*H[j]*detJ*dms;
                            kuc += ((mat3dd(1.0)*phif - D[jsol]/d0[jsol])*BF.force(mp)*mt.m_dkdc[jsol][isol]*mt.m_c[jsol]*M[jsol]*H[i]*H[j] + (D[jsol]*d0p[jsol][isol]/d0[jsol]/d0[jsol] - dDdc[jsol][isol]/d0[jsol])*BF.force(mp)*M[jsol]*mt.m_k[jsol]*mt.m_c[jsol]*H[i]*H[j])*dms*detJ;
                            kcc = -(D[isol]*mt.m_dkdc[isol][jsol]*mt.m_c[isol] + dDdc[isol][jsol]*mt.m_k[isol]*mt.m_c[isol])*BF.force(mp)*gradN[i]*M[isol]/R/T*H[j]*phif*detJ;
                            for (int ksol=0; ksol<nsol; ++ksol)
                            {
                                if(jsol==ksol)
                                {
                                    kcc += -(D[jsol]*mt.m_dkdc[jsol][jsol]*mt.m_c[jsol] + D[jsol]*mt.m_k[jsol] + dDdc[jsol][jsol]*mt.m_k[jsol]*mt.m_c[jsol])*BF.force(mp)*gradN[i]*M[jsol]*z[jsol]*penalty/R/T*H[j]*phif*detJ;
                                }
                                else
                                {
                                    kcc += -(D[ksol]*mt.m_dkdc[ksol][jsol]*mt.m_c[ksol] + dDdc[ksol][jsol]*mt.m_k[ksol]*mt.m_c[ksol])*BF.force(mp)*gradN[i]*M[ksol]*z[ksol]*penalty/R/T*H[j]*phif*detJ;
                                }
                            }
                        }
                        ke[i7+7+isol][j7+7+jsol] += kcc;
                    }
                    ke[i7+0][j7+7+isol] += kuc.x;
                    ke[i7+1][j7+7+isol] += kuc.y;
                    ke[i7+2][j7+7+isol] += kuc.z;
                    ke[i7+3][j7+7+isol] += kwc.x;
                    ke[i7+4][j7+7+isol] += kwc.y;
                    ke[i7+5][j7+7+isol] += kwc.z;
                    ke[i7+7+isol][j7+0] += kcu.x;
                    ke[i7+7+isol][j7+1] += kcu.y;
                    ke[i7+7+isol][j7+2] += kcu.z;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEMultiphasicFSIDomain3D::ElementStiffness(FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i7, j, j7, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    const int nsol = m_pMat->Solutes();
    const int ndpn = 7 + nsol;
    
    const int nreact = m_pMat->Reactions();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    vector<mat3d> gradgradN(neln);
    
    double dt = tp.timeIncrement;
    double a = tp.gamma/(tp.beta*dt);
    double c = tp.alpham/(tp.alphaf*tp.gamma*dt);
    
    double dtrans = m_btrans ? 1 : m_sseps;
    
    double *H, *Gr, *Gs, *Gt;
    vec3d g[3], dg[3][3];
    
    // jacobian
    double Ji[3][3], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjact(el, Ji, n, tp.alphaf)*gw[n]*tp.alphaf;
        
        ContraBaseVectors(el, n, g, tp.alphaf);
        ContraBaseVectorDerivatives(el, n, dg, tp.alphaf);
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // get the shape function derivatives
        double* Grr = el.Grr(n); double* Grs = el.Grs(n); double* Grt = el.Grt(n);
        double* Gsr = el.Gsr(n); double* Gss = el.Gss(n); double* Gst = el.Gst(n);
        double* Gtr = el.Gtr(n); double* Gts = el.Gts(n); double* Gtt = el.Gtt(n);
        
        // setup the material point
        // NOTE: deformation gradient and determinant have already been evaluated in the stress routine
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& et = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFSIMaterialPoint& fpt = *(mp.ExtractData<FEFSIMaterialPoint>());
        FEBiphasicFSIMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
        FEMultiphasicFSIMaterialPoint& mt = *(mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
        double Jf = 1 + pt.m_ef;
        
        // get the tangents
        mat3ds se = m_pMat->Solid()->Stress(mp);
        tens4ds cs = m_pMat->Solid()->Tangent(mp);
        mat3ds sv = m_pMat->Fluid()->GetViscous()->Stress(mp);
        mat3ds svJ = m_pMat->Fluid()->GetViscous()->Tangent_Strain(mp);
        tens4ds cv = m_pMat->Fluid()->Tangent_RateOfDeformation(mp);
        double pa = m_pMat->PressureActual(mp);
        double dp = m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        double d2p = m_pMat->Fluid()->Tangent_Pressure_Strain_Strain(mp);
        vec3d gradp = pt.m_gradef*dp;
        // Jsdot/Js = div(vs)
        double dJsoJ = fpt.m_Jdot/et.m_J;
        mat3ds km1 = m_pMat->InvPermeability(mp);
        
        //Include dependence of permeability on displacement
        tens4dmm K = m_pMat->Permeability_Tangent(mp);
        
        double phif = m_pMat->Porosity(mp);
        double phis = m_pMat->SolidVolumeFrac(mp);
        vec3d gradphif = m_pMat->gradPorosity(mp);
        vec3d gradphifphis = m_pMat->gradPhifPhis(mp);
        
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double dms = m_pMat->m_diffMtmSupp;
        double penalty = m_pMat->m_penalty;
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        double dodJ = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Strain(mp);
        vector<double> MM(nsol);
        vector<int> z(nsol);
        vector<double> d0(nsol);
        vector<vector<double>> d0p(nsol, vector<double>(nsol));
        vector<mat3ds> D(nsol);
        vector<mat3ds> Dm1(nsol);
        vector<tens4dmm> dDdE(nsol);
        vector<vector<mat3ds>> dDdc(nsol, vector<mat3ds>(nsol));
        vector<tens4d> Dm1Dm1(nsol);
        vector<tens4d> dDdEfull(nsol);
        vector<vec3d> flux(nsol);
        vector<double> dodc(nsol);
        
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            MM[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            D[isol] = m_pMat->Diffusivity(mp, isol);
            Dm1[isol] = m_pMat->InvDiffusivity(mp, isol);
            Dm1Dm1[isol] = dyad2(Dm1[isol],Dm1[isol]);
            dDdE[isol] = m_pMat->Diffusivity_Tangent_Strain(mp, isol);
            dDdEfull[isol] = tens4d(dDdE[isol]);
            flux[isol] = -mt.m_gradc[isol]*phif + fpt.m_w*mt.m_c[isol]/d0[isol];
            dodc[isol] = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp, isol);
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                d0p[isol][jsol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, jsol);
                dDdc[isol][jsol] = m_pMat->Diffusivity_Tangent_Concentration(mp, isol, jsol);
            }
        }
        
        vector<double> dkdt(nsol,0);
        for (int isol=0; isol<nsol; ++isol)
        {
            dkdt[isol] = mt.m_dkdJ[isol]*fpt.m_Jdot;
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                dkdt[isol] += mt.m_dkdc[isol][jsol]*mt.m_cdot[jsol];
            }
        }
        
        // evaluate the chat
        double vbarzeta = 0.0;
        vector<double> vzeta(nsol,0.0);
        mat3ds vbardzdE = mat3ds(0.0);
        vector<mat3ds> vdzdE(nsol,mat3ds(0.0));
        vector<double> vbardzdc(nsol,0.0);
        vector<vector<double>> vdzdc(nsol,vector<double>(nsol, 0.0));
        
        // chemical reactions
        for (i=0; i<nreact; ++i) {
            double vbar = m_pMat->GetReaction(i)->m_Vbar;
            mat3ds dzdE = m_pMat->GetReaction(i)->Tangent_ReactionSupply_Strain(mp);
            double zeta = m_pMat->GetReaction(i)->ReactionSupply(mp);
            vbarzeta += vbar*zeta;
            vbardzdE += dzdE*vbar;
            for (int isol = 0; isol < nsol; ++isol)
            {
                double v = m_pMat->GetReaction(i)->m_v[isol];
                double dzdc1 = m_pMat->GetReaction(i)->Tangent_ReactionSupply_Concentration(mp,isol);
                vzeta[isol] += v*zeta;
                vdzdE[isol] += dzdE*v;
                vbardzdc[isol] += vbar*dzdc1;
                for (int jsol = 0; jsol<nsol; ++jsol)
                {
                    double dzdc2 = m_pMat->GetReaction(i)->Tangent_ReactionSupply_Concentration(mp,jsol);
                    vdzdc[isol][jsol] += v*dzdc2;
                }
            }
        }
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate spatial gradgrad of shape functions
        for (i=0; i<neln; ++i) {
            gradgradN[i] = (((dg[0][0] & g[0]) + (dg[0][1] & g[1]) + (dg[0][2] & g[2]))*Gr[i]
                            + ((dg[1][0] & g[0]) + (dg[1][1] & g[1]) + (dg[1][2] & g[2]))*Gs[i]
                            + ((dg[2][0] & g[0]) + (dg[2][1] & g[1]) + (dg[2][2] & g[2]))*Gt[i]
                            + (g[0] & g[0])*Grr[i] + (g[0] & g[1])*Gsr[i] + (g[0] & g[2])*Gtr[i]
                            + (g[1] & g[0])*Grs[i] + (g[1] & g[1])*Gss[i] + (g[1] & g[2] )*Gts[i]
                            + (g[2] & g[0])*Grt[i] + (g[2] & g[1])*Gst[i] + (g[2] & g[2])*Gtt[i]);
        }
        
        // evaluate stiffness matrix
        for (i=0, i7=0; i<neln; ++i, i7 += ndpn)
        {
            for (j=0, j7 = 0; j<neln; ++j, j7 += ndpn)
            {
                mat3d M = mat3dd(a*dtrans) - et.m_L;
                tens4d km1km1 = dyad2(km1,km1);
                tens4d Kfull = tens4d(K);
                
                mat3d Kuu = (sv*((gradN[i]&gradN[j])*phis/phif + (gradN[j]&gradN[i]))*phis - vdotTdotv(gradN[i], cv, gradN[j])*(-bpt.m_Lw.sym()*phis/(phif*phif) + M)*phis - vdotTdotv(gradN[i], cv, fpt.m_w) * ((gradphif&gradN[j])*2.0*phis/phif + (gradN[j]&gradphif))*phis/(phif*phif) + vdotTdotv(gradN[i], cv, fpt.m_w) * ((-gradphif&gradN[j]) + gradgradN[j]*phis)*phis/(phif*phif) + mat3dd((se*gradN[i])*gradN[j]) + vdotTdotv(gradN[i], cs, gradN[j]) -  sv*((gradphifphis&gradN[j])*2.0*phis*phis/phif + (gradN[j]&gradphifphis)*phis - gradgradN[j])*H[i]*phis/phif + vdotTdotv(gradphifphis, cv, gradN[j])*(-bpt.m_Lw.sym()*phis/(phif*phif) + M)*phis*phis/phif*H[i] + vdotTdotv(gradphifphis, cv, fpt.m_w) * ((gradphif&gradN[j])*2.0*phis/phif + (gradN[j]&gradphif))*phis*phis/(phif*phif*phif)*H[i] - vdotTdotv(gradphifphis, cv, fpt.m_w) * ((-gradphif&gradN[j]) + gradgradN[j]*phis)*phis*phis/(phif*phif*phif)*H[i] + (-((km1*fpt.m_w)&gradN[j])*2.0 + (gradN[j]&(km1*fpt.m_w)) + km1*(gradN[j]*fpt.m_w) + ddot(ddot(km1km1,Kfull),mat3dd(gradN[j]*fpt.m_w)))*H[i])*detJ; //Old Nat BC
                
                //mat3d Kuu = (sv*((gradN[i]&gradN[j])*phis/phif + (gradN[j]&gradN[i]))*phis - vdotTdotv(gradN[i], cv, gradN[j])*(-bpt.m_Lw.sym()*phis/(phif*phif) + M)*phis - vdotTdotv(gradN[i], cv, fpt.m_w) * ((gradphif&gradN[j])*2.0*phis/phif + (gradN[j]&gradphif))*phis/(phif*phif) + vdotTdotv(gradN[i], cv, fpt.m_w) * ((-gradphif&gradN[j]) + gradgradN[j]*phis)*phis/(phif*phif) + mat3dd((se*gradN[i])*gradN[j]) + vdotTdotv(gradN[i], cs, gradN[j]) -  sv*((gradphifphis&gradN[j])*2.0*phis*phis/phif + (gradN[j]&gradphifphis)*phis - gradgradN[j])*H[i]*phis/phif + vdotTdotv(gradphifphis, cv, gradN[j])*(-bpt.m_Lw.sym()*phis/(phif*phif) + M)*phis*phis/phif*H[i] + vdotTdotv(gradphifphis, cv, fpt.m_w) * ((gradphif&gradN[j])*2.0*phis/phif + (gradN[j]&gradphif))*phis*phis/(phif*phif*phif)*H[i] - vdotTdotv(gradphifphis, cv, fpt.m_w) * ((-gradphif&gradN[j]) + gradgradN[j]*phis)*phis*phis/(phif*phif*phif)*H[i] + (-((km1*fpt.m_w)&gradN[j])*2.0 + (gradN[j]&(km1*fpt.m_w)) + km1*(gradN[j]*fpt.m_w) + ddot(ddot(km1km1,Kfull),mat3dd(gradN[j]*fpt.m_w)))*H[i] - ((gradp&gradN[j])-(gradN[j]&gradp))*H[i] - (gradN[i]&gradN[j])*pa + (gradN[j]&gradN[i])*pa)*detJ; //New Nat BC
                
                mat3d Kuw = (vdotTdotv(gradN[i], cv, (gradphif*H[j]/phif-gradN[j]))*phis/phif + vdotTdotv(gradphifphis, cv, (-gradphif*H[j]/phif + gradN[j]))*H[i]*phis*phis/(phif*phif) - km1*H[i]*H[j])*detJ;
                
                vec3d kuJ = ((-svJ*gradN[i])*H[j]*phis + svJ*gradphifphis*H[j]*H[i]*phis*phis/phif)*detJ; //Old Nat BC
                //vec3d kuJ = (-(mat3dd(1.0)*dp+svJ*phis)*gradN[i]*H[j] + ((-pt.m_gradef*d2p + svJ*gradphifphis*phis*phis/phif)*H[j] - gradN[j]*dp)*H[i])*detJ; //New Nat BC
                
                mat3d Kwu = (((gradp&gradN[j])-(gradN[j]&gradp))*H[i] + sv*((gradphif&gradN[j])*(2*phis/phif)+(gradN[j]&gradphif) - gradgradN[j]*phis)*(H[i]/phif) - vdotTdotv(gradphif, cv, gradN[j])*(-bpt.m_Lw.sym()*phis/(phif*phif) + M)*H[i]/phif - vdotTdotv(gradphif, cv, fpt.m_w) * ((gradphif&gradN[j])*2.0*phis/phif + (gradN[j]&gradphif))*H[i]/(phif*phif*phif) + vdotTdotv(gradphif, cv, fpt.m_w) * (-(gradphif&gradN[j]) + gradgradN[j]*phis)*H[i]/(phif*phif*phif) + sv*((gradN[i]&gradN[j])*(1-phis/phif)-(gradN[j]&gradN[i])) + vdotTdotv(gradN[i], cv, gradN[j])*(-bpt.m_Lw.sym()*phis/(phif*phif) + M) + vdotTdotv(gradN[i], cv, fpt.m_w) * ((gradphif&gradN[j])*2.0*phis/phif + (gradN[j]&gradphif))/(phif*phif) + vdotTdotv(gradN[i], cv, fpt.m_w) * ((gradphif&gradN[j]) - gradgradN[j]*phis)/(phif*phif) + (((km1*fpt.m_w)&gradN[j])*2.0 - (gradN[j]&(km1*fpt.m_w)) - km1*(gradN[j]*fpt.m_w) - ddot(ddot(km1km1,Kfull),mat3dd(gradN[j]*fpt.m_w)))*H[i])*detJ;
                
                mat3d Kww = ((vdotTdotv(gradphif, cv, (gradphif*H[j]/phif-gradN[j]))/(phif*phif) + km1*H[j])*H[i] + vdotTdotv(gradN[i], cv, (-gradphif*H[j]/phif+gradN[j]))/phif)*detJ;
                
                vec3d kwJ = ((svJ*gradN[i])*H[j] +(gradN[j]*dp+(pt.m_gradef*d2p - svJ*gradphif/phif)*H[j])*H[i])*detJ;
                
                vec3d kJu = (((gradN[j]&fpt.m_w) - mat3dd(gradN[j]*fpt.m_w)) * gradN[i] + ((gradN[j]*pt.m_efdot + ((gradN[j]&fpt.m_w) - mat3dd(gradN[j]*fpt.m_w))*pt.m_gradef)/Jf - gradN[j]*(dJsoJ + a*dtrans) + et.m_L.transpose()*gradN[j]*dtrans)*H[i] + (mat3dd(1.0)*vbarzeta + et.m_F*vbardzdE*et.m_F.transpose()*phif)*gradN[j]*H[i])*detJ;
                
                vec3d kJw = ((pt.m_gradef*(H[i]/Jf) + gradN[i])*H[j])*detJ;
                
                double kJJ = ((c*phif*dtrans - (pt.m_efdot*phif + pt.m_gradef*fpt.m_w)/Jf)*H[j] + gradN[j]*fpt.m_w)*H[i]/Jf*detJ;
                
                for (int isol = 0; isol < nsol; ++isol)
                {
                    Kuu += (-(gradN[i]&gradN[j])*(osmc + et.m_J*dodJ)*R*T*mt.m_k[isol]*mt.m_c[isol] + (-(gradN[i]&gradN[j])*et.m_J*mt.m_dkdJ[isol] + ((gradN[j]&gradN[i]))*mt.m_k[isol])*mt.m_c[isol]*R*T*osmc - (mt.m_gradc[isol]&gradN[j])*mt.m_k[isol]*R*T*H[i] + (-(mt.m_gradc[isol]&gradN[j])*et.m_J*mt.m_dkdJ[isol] + (gradN[j]&mt.m_gradc[isol])*mt.m_k[isol])*H[i]*phif*R*T + (-((Dm1[isol]*mt.m_j[isol])&gradN[j])*2.0 + (gradN[j]&(Dm1[isol]*mt.m_j[isol])) + Dm1[isol]*(gradN[j]*mt.m_j[isol]) + ddot(ddot(Dm1Dm1[isol],dDdEfull[isol]),mat3dd(gradN[j]*mt.m_j[isol])))*H[i]*R*T + Dm1[isol]*(mat3dd((D[isol]*flux[isol])*gradN[j]) - D[isol]*(flux[isol]&gradN[j]) - dDdEfull[isol].dot(mat3dd(gradN[j]*flux[isol])) - D[isol]*(gradN[j]&flux[isol]))*mt.m_k[isol]*H[i]*R*T - (flux[isol]&gradN[j])*et.m_J*mt.m_dkdJ[isol]*H[i]*R*T - (-(mt.m_gradc[isol]&gradN[j])*phis + (gradN[j]&mt.m_gradc[isol])*phif)*mt.m_k[isol]*H[i]*R*T + (mt.m_j[isol]&gradN[j])*H[i]*(1/phif - phis/(phif*phif))*R*T/d0[isol] + (mat3dd((D[isol]*flux[isol])*gradN[j]) - D[isol]*(flux[isol]&gradN[j]) + dDdEfull[isol].dot(mat3dd(gradN[j]*flux[isol])) + D[isol]*(gradN[j]&flux[isol]))*H[i]*R*T/phif/d0[isol]*mt.m_k[isol] + D[isol]*(flux[isol]&gradN[j])*H[i]*R*T/phif*et.m_J*mt.m_dkdJ[isol]/d0[isol] + D[isol]*(-(mt.m_gradc[isol]&gradN[j])*phis + (gradN[j]&mt.m_gradc[isol])*phif)*H[i]*R*T/phif*mt.m_k[isol]/d0[isol] + (fpt.m_w&gradN[j])*H[i]*(phis/phif*mt.m_k[isol] - et.m_J*mt.m_dkdJ[isol])*phis/phif*R*T*mt.m_c[isol]/d0[isol])*dms*detJ; //Old Nat BC
                    //Kuu += (-(gradN[i]&gradN[j])*et.m_J*dodJ*R*T*mt.m_k[isol]*mt.m_c[isol] -(gradN[i]&gradN[j])*et.m_J*mt.m_dkdJ[isol]*mt.m_c[isol]*R*T*osmc - (mt.m_gradc[isol]&gradN[j])*mt.m_k[isol]*R*T*H[i] + (-(mt.m_gradc[isol]&gradN[j])*et.m_J*mt.m_dkdJ[isol] + (gradN[j]&mt.m_gradc[isol])*mt.m_k[isol])*H[i]*phif*R*T + (-((Dm1[isol]*mt.m_j[isol])&gradN[j])*2.0 + (gradN[j]&(Dm1[isol]*mt.m_j[isol])) + Dm1[isol]*(gradN[j]*mt.m_j[isol]) + ddot(ddot(Dm1Dm1[isol],dDdEfull[isol]),mat3dd(gradN[j]*mt.m_j[isol])))*H[i]*R*T + Dm1[isol]*(mat3dd((D[isol]*flux[isol])*gradN[j]) - D[isol]*(flux[isol]&gradN[j]) - dDdEfull[isol].dot(mat3dd(gradN[j]*flux[isol])) - D[isol]*(gradN[j]&flux[isol]))*mt.m_k[isol]*H[i]*R*T - (flux[isol]&gradN[j])*et.m_J*mt.m_dkdJ[isol]*H[i]*R*T - (-(mt.m_gradc[isol]&gradN[j])*phis + (gradN[j]&mt.m_gradc[isol])*phif)*mt.m_k[isol]*H[i]*R*T + (mt.m_j[isol]&gradN[j])*H[i]*(1/phif - phis/(phif*phif))*R*T/d0[isol] + (mat3dd((D[isol]*flux[isol])*gradN[j]) - D[isol]*(flux[isol]&gradN[j]) + dDdEfull[isol].dot(mat3dd(gradN[j]*flux[isol])) + D[isol]*(gradN[j]&flux[isol]))*H[i]*R*T/phif/d0[isol]*mt.m_k[isol] + D[isol]*(flux[isol]&gradN[j])*H[i]*R*T/phif*et.m_J*mt.m_dkdJ[isol]/d0[isol] + D[isol]*(-(mt.m_gradc[isol]&gradN[j])*phis + (gradN[j]&mt.m_gradc[isol])*phif)*H[i]*R*T/phif*mt.m_k[isol]/d0[isol] + (fpt.m_w&gradN[j])*H[i]*(phis/phif*mt.m_k[isol] - et.m_J*mt.m_dkdJ[isol])*phis/phif*R*T*mt.m_c[isol]/d0[isol])*dms*detJ; //New Nat BC
                    
                    Kwu += ((fpt.m_w&gradN[j])*H[i]*(1.0/phif-phis/(phif*phif))*R*T*mt.m_k[isol]*mt.m_c[isol]/d0[isol] + (fpt.m_w&gradN[j])*H[i]/phif*et.m_J*R*T*mt.m_dkdJ[isol]*mt.m_c[isol]/d0[isol] - (mt.m_j[isol]&gradN[j])*H[i]*(1.0/phif-phis/(phif*phif))*R*T/d0[isol] - (mat3dd((D[isol]*flux[isol])*gradN[j]) - D[isol]*(flux[isol]&gradN[j]) + dDdEfull[isol].dot(mat3dd(flux[isol]*gradN[j])) + D[isol]*(gradN[j]&flux[isol]) + D[isol]*(-(mt.m_gradc[isol]&gradN[j])*phis + (gradN[j]&mt.m_gradc[isol])*phif))*H[i]*R*T/phif/d0[isol]*mt.m_k[isol] - D[isol]*(flux[isol]&gradN[j])*et.m_J*mt.m_dkdJ[isol]*H[i]*R*T/phif/d0[isol])*detJ*dms;
                    
                    Kww += (mat3dd(1.0) - D[isol]/d0[isol])*mt.m_k[isol]*mt.m_c[isol]/d0[isol]*R*T/phif*H[i]*H[j]*dms*detJ;
                    
                    Kuw += -((mat3dd(1.0) - D[isol]/d0[isol]/phif)*mt.m_k[isol]*mt.m_c[isol]/d0[isol]*R*T*H[i]*H[j] + mat3dd(1.0)*mt.m_k[isol]*mt.m_c[isol]/d0[isol]*R*T*phis/phif*H[i]*H[j])*dms*detJ;
                }
                
                ke[i7  ][j7  ] += Kuu(0,0); ke[i7  ][j7+1] += Kuu(0,1); ke[i7  ][j7+2] += Kuu(0,2);
                ke[i7+1][j7  ] += Kuu(1,0); ke[i7+1][j7+1] += Kuu(1,1); ke[i7+1][j7+2] += Kuu(1,2);
                ke[i7+2][j7  ] += Kuu(2,0); ke[i7+2][j7+1] += Kuu(2,1); ke[i7+2][j7+2] += Kuu(2,2);
                
                ke[i7  ][j7+3] += Kuw(0,0); ke[i7  ][j7+4] += Kuw(0,1); ke[i7  ][j7+5] += Kuw(0,2);
                ke[i7+1][j7+3] += Kuw(1,0); ke[i7+1][j7+4] += Kuw(1,1); ke[i7+1][j7+5] += Kuw(1,2);
                ke[i7+2][j7+3] += Kuw(2,0); ke[i7+2][j7+4] += Kuw(2,1); ke[i7+2][j7+5] += Kuw(2,2);
                
                ke[i7+3][j7  ] += Kwu(0,0); ke[i7+3][j7+1] += Kwu(0,1); ke[i7+3][j7+2] += Kwu(0,2);
                ke[i7+4][j7  ] += Kwu(1,0); ke[i7+4][j7+1] += Kwu(1,1); ke[i7+4][j7+2] += Kwu(1,2);
                ke[i7+5][j7  ] += Kwu(2,0); ke[i7+5][j7+1] += Kwu(2,1); ke[i7+5][j7+2] += Kwu(2,2);
                
                ke[i7+3][j7+3] += Kww(0,0); ke[i7+3][j7+4] += Kww(0,1); ke[i7+3][j7+5] += Kww(0,2);
                ke[i7+4][j7+3] += Kww(1,0); ke[i7+4][j7+4] += Kww(1,1); ke[i7+4][j7+5] += Kww(1,2);
                ke[i7+5][j7+3] += Kww(2,0); ke[i7+5][j7+4] += Kww(2,1); ke[i7+5][j7+5] += Kww(2,2);
                
                ke[i7+0][j7+6] += kuJ.x;
                ke[i7+1][j7+6] += kuJ.y;
                ke[i7+2][j7+6] += kuJ.z;
                
                ke[i7+3][j7+6] += kwJ.x;
                ke[i7+4][j7+6] += kwJ.y;
                ke[i7+5][j7+6] += kwJ.z;
                
                ke[i7+6][j7  ] += kJu.x;
                ke[i7+6][j7+1] += kJu.y;
                ke[i7+6][j7+2] += kJu.z;
                
                ke[i7+6][j7+3] += kJw.x;
                ke[i7+6][j7+4] += kJw.y;
                ke[i7+6][j7+5] += kJw.z;
                
                ke[i7+6][j7+6] += kJJ;
                
                for (int isol=0; isol<nsol; ++isol)
                {
                    vec3d kcw = D[isol]*gradN[i]/d0[isol]*mt.m_k[isol]*mt.m_c[isol]*H[j]*detJ;
                    
                    vec3d kcu = (((gradN[j]&mt.m_j[isol]) - mat3dd(mt.m_j[isol]*gradN[j]))*gradN[i] + (mat3dd((D[isol]*flux[isol])*gradN[j]) - (gradN[j]&(D[isol]*flux[isol])) + (dDdEfull[isol].dot(mat3dd(gradN[j]*flux[isol]))).transpose() + (flux[isol]&gradN[j])*D[isol])*gradN[i]*mt.m_k[isol] + gradN[j]*((D[isol]*flux[isol])*gradN[i])*et.m_J*mt.m_dkdJ[isol] + (-(gradN[j]&mt.m_gradc[isol])*phis + (mt.m_gradc[isol]&gradN[j])*phif)*D[isol]*gradN[i]*mt.m_k[isol] + (mat3dd(vzeta[isol]) + et.m_F*vdzdE[isol]*et.m_F.transpose()*phif)*gradN[j]*H[i] - (mat3dd(2*mt.m_dkdJ[isol]*fpt.m_Jdot*mt.m_c[isol]) + (mat3dd(dJsoJ + a*dtrans) - et.m_L.transpose())*mt.m_k[isol]*mt.m_c[isol])*gradN[j]*H[i] - (mat3dd(dJsoJ + a*dtrans) - et.m_L.transpose())*gradN[j]*H[i]*et.m_J*phif*mt.m_c[isol]*mt.m_dkdJ[isol] - gradN[j]*H[i]*(mt.m_k[isol] + et.m_J*phif*mt.m_dkdJ[isol])*mt.m_cdot[isol])*detJ;
                    
                    vec3d kwc = vec3d(0);
                    
                    vec3d kuc = vec3d(0);
                    
                    double kJc = H[i]*H[j]*phif*vbardzdc[isol]*mt.m_k[isol]*detJ;
                    
                    for(int jsol=0; jsol<nsol; ++jsol)
                    {
                        kcw += D[jsol]*gradN[i]/d0[jsol]*mt.m_k[jsol]*mt.m_c[jsol]*z[jsol]*H[j]*detJ;
                        
                        kcu += (((gradN[j]&mt.m_j[jsol])*z[jsol]*penalty - mat3dd(mt.m_j[isol]*gradN[j]*z[jsol]*penalty))*gradN[i] + (mat3dd((D[jsol]*flux[jsol])*gradN[j]) - (gradN[j]&(D[jsol]*flux[jsol])) + (dDdEfull[jsol].dot(mat3dd(gradN[j]*flux[jsol]))).transpose() + (flux[jsol]&gradN[j])*D[jsol])*gradN[i]*mt.m_k[jsol]*z[jsol]*penalty + gradN[j]*((D[jsol]*flux[jsol])*gradN[i])*et.m_J*mt.m_dkdJ[jsol]*z[jsol]*penalty + (-(gradN[j]&mt.m_gradc[jsol])*phis + (mt.m_gradc[jsol]&gradN[j])*phif)*D[jsol]*gradN[i]*mt.m_k[jsol]*z[jsol]*penalty + mat3dd(vdzdc[isol][jsol]*mt.m_dkdJ[jsol]*mt.m_c[jsol])*gradN[j]*phif*et.m_J*H[i] - gradN[j]*H[i]*mt.m_dkdc[isol][jsol]*mt.m_c[isol]*mt.m_cdot[jsol])*detJ;
                        
                        kJc += H[i]*H[j]*phif*vbardzdc[isol]*mt.m_dkdc[jsol][isol]*mt.m_c[jsol]*detJ;
                        
                        double kcc = 0;
                        
                        if (isol == jsol)
                        {
                            kwc += (fpt.m_w*H[i]*H[j]*R*T/phif*(mt.m_dkdc[isol][isol]*mt.m_c[isol]/d0[isol] + mt.m_k[isol]/d0[isol] - mt.m_k[isol]*mt.m_c[isol]*d0p[isol][isol]/(d0[isol]*d0[isol])) + mt.m_j[isol]*H[i]*H[j]*R*T/phif*d0p[isol][isol]/(d0[isol]*d0[isol]) - (dDdc[isol][isol]*mt.m_k[isol] + D[isol]*mt.m_dkdc[isol][isol])*flux[isol]*H[i]*H[j]*R*T/phif/d0[isol] - D[isol]*(-gradN[j]*phif + fpt.m_w*(1.0/d0[isol] - mt.m_c[isol]*d0p[isol][isol]/(d0[isol]*d0[isol]))*H[j])*H[i]*R*T/phif/d0[isol]*mt.m_k[isol])*detJ*dms;
                            
                            kuc += (-gradN[i]*H[j]*R*T*(dodc[isol]*mt.m_k[isol]*mt.m_c[isol] + osmc*mt.m_dkdc[isol][isol]*mt.m_c[isol] + osmc*mt.m_k[isol]) - (mt.m_gradc[isol]*mt.m_dkdc[isol][isol]*H[j] + gradN[j]*mt.m_k[isol])*H[i]*R*T*phif + (Dm1Dm1[isol].dot(dDdc[isol][isol]) - mat3dd(d0p[isol][isol]/phif/d0[isol]/d0[isol]))*mt.m_j[isol]*H[i]*H[j]*R*T - (Dm1[isol] - mat3dd(1.0/phif/d0[isol]))*(dDdc[isol][isol]*mt.m_k[isol] + D[isol]*mt.m_dkdc[isol][isol])*flux[isol]*H[i]*H[j]*R*T - (mat3dd(1.0) - D[isol]/d0[isol]/phif)*(-gradN[j]*phif + fpt.m_w*H[j]*(1.0/d0[isol] - mt.m_c[isol]*d0p[isol][isol]/(d0[isol]*d0[isol])))*H[i]*R*T*mt.m_k[isol] - fpt.m_w*H[i]*H[j]*R*T*phis/phif*(mt.m_dkdc[isol][isol]*mt.m_c[isol]/d0[isol] + mt.m_k[isol]/d0[isol] - mt.m_k[isol]*mt.m_c[isol]*d0p[isol][isol]/(d0[isol]*d0[isol])))*detJ*dms;
                            
                            kcc += (((dDdc[isol][isol]*mt.m_k[isol] + D[isol]*mt.m_dkdc[isol][isol])*flux[isol])*gradN[i]*H[j] + (D[isol]*(-gradN[j]*phif + fpt.m_w*H[j]*(1.0/d0[isol] - mt.m_c[isol]*d0p[isol][isol]/(d0[isol]*d0[isol]))))*gradN[i]*mt.m_k[isol] + H[i]*H[j]*phif*vdzdc[isol][isol]*mt.m_k[isol] - H[i]*H[j]*dJsoJ*(mt.m_dkdc[isol][isol]*mt.m_c[isol] + mt.m_k[isol]) - H[i]*H[j]*phif*mt.m_dkdJ[isol]*fpt.m_Jdot - H[i]*H[j]*phif*mt.m_c[isol]*mt.m_dkdc[isol][isol]*c*dtrans - H[i]*H[j]*phif*(mt.m_dkdc[isol][isol]*mt.m_cdot[isol] + mt.m_k[isol]*c*dtrans))*detJ;
                            
                            for (int ksol=0; ksol<nsol; ++ksol)
                            {
                                if(isol==ksol)
                                {
                                    kcc += (((dDdc[isol][isol]*mt.m_k[isol] + D[isol]*mt.m_dkdc[isol][isol])*flux[isol])*gradN[i]*H[j]*z[isol]*penalty + (D[isol]*(-gradN[j]*phif + fpt.m_w*H[j]*(1.0/d0[isol] - mt.m_c[isol]*d0p[isol][isol]/(d0[isol]*d0[isol]))))*gradN[i]*mt.m_k[isol]*z[isol]*penalty + H[i]*H[j]*phif*vdzdc[isol][isol]*mt.m_dkdc[isol][isol]*mt.m_c[isol] - H[i]*H[j]*phif*mt.m_dkdc[isol][isol]*mt.m_cdot[isol])*detJ;
                                }
                                else
                                {
                                    kcc += (((dDdc[ksol][isol]*mt.m_k[ksol] + D[ksol]*mt.m_dkdc[ksol][isol])*flux[ksol])*gradN[i]*H[j]*z[ksol]*penalty - (D[ksol]*fpt.m_w*H[j]*mt.m_c[ksol]*d0p[ksol][isol]/(d0[ksol]*d0[ksol]))*gradN[i]*mt.m_k[ksol]*z[ksol]*penalty + H[i]*H[j]*phif*vdzdc[isol][isol]*mt.m_dkdc[ksol][isol]*mt.m_c[ksol] - H[i]*H[j]*phif*mt.m_dkdc[isol][ksol]*mt.m_cdot[ksol])*detJ;
                                }
                            }
                        }
                        else
                        {
                            kwc += (fpt.m_w*H[i]*H[j]*R*T/phif*(mt.m_dkdc[jsol][isol]*mt.m_c[jsol]/d0[jsol] - mt.m_k[jsol]*mt.m_c[jsol]*d0p[jsol][isol]/(d0[jsol]*d0[jsol])) + mt.m_j[jsol]*H[i]*H[j]*R*T/phif*d0p[jsol][isol]/(d0[jsol]*d0[jsol]) - (dDdc[jsol][isol]*mt.m_k[jsol] + D[jsol]*mt.m_dkdc[jsol][isol])*flux[jsol]*H[i]*H[j]*R*T/phif/d0[jsol] + D[jsol]*fpt.m_w*mt.m_k[jsol]*mt.m_c[jsol]*d0p[jsol][isol]/(d0[jsol]*d0[jsol]*d0[jsol])*H[j]*H[i]*R*T/phif)*detJ*dms;
                            
                            kuc += (-gradN[i]*H[j]*R*T*(dodc[isol]*mt.m_k[jsol]*mt.m_c[jsol] + osmc*mt.m_dkdc[jsol][isol]*mt.m_c[jsol]) - mt.m_gradc[jsol]*mt.m_dkdc[jsol][isol]*H[i]*H[j]*R*T*phif + (Dm1Dm1[jsol].dot(dDdc[jsol][isol]) - mat3dd(d0p[jsol][isol]/phif/d0[jsol]/d0[jsol]))*mt.m_j[jsol]*H[i]*H[j]*R*T - (Dm1[jsol] - mat3dd(1.0/phif/d0[jsol]))*(dDdc[jsol][isol]*mt.m_k[jsol] + D[jsol]*mt.m_dkdc[jsol][isol])*flux[jsol]*H[i]*H[j]*R*T + (mat3dd(1.0) - D[jsol]/d0[jsol]/phif)*fpt.m_w*H[j]*mt.m_c[jsol]*d0p[jsol][isol]/(d0[jsol]*d0[jsol])*H[i]*R*T*mt.m_k[jsol] - fpt.m_w*H[i]*H[j]*R*T*phis/phif*(mt.m_dkdc[jsol][isol]*mt.m_c[jsol]/d0[jsol] - mt.m_k[jsol]*mt.m_c[jsol]*d0p[jsol][isol]/(d0[jsol]*d0[jsol])))*detJ*dms;
                            
                            kcc += (((dDdc[isol][jsol]*mt.m_k[isol] + D[isol]*mt.m_dkdc[isol][jsol])*flux[isol])*gradN[i]*H[j] - (D[isol]*fpt.m_w*H[j]*mt.m_c[isol]*d0p[isol][jsol]/(d0[isol]*d0[isol]))*gradN[i]*mt.m_k[isol] + H[i]*H[j]*phif*vdzdc[isol][jsol]*mt.m_k[jsol] - H[i]*H[j]*dJsoJ*mt.m_dkdc[isol][jsol]*mt.m_c[isol] - H[i]*H[j]*phif*mt.m_c[isol]*mt.m_dkdc[isol][jsol]*c*dtrans - H[i]*H[j]*phif*mt.m_dkdc[isol][jsol]*mt.m_cdot[isol])*detJ;
                            
                            for (int ksol=0; ksol<nsol; ++ksol)
                            {
                                if(jsol==ksol)
                                {
                                    kcc += (((dDdc[jsol][jsol]*mt.m_k[jsol] + D[jsol]*mt.m_dkdc[jsol][jsol])*flux[jsol])*gradN[i]*H[j]*z[jsol]*penalty + (D[jsol]*(-gradN[j]*phif + fpt.m_w*H[j]*(1.0/d0[jsol] - mt.m_c[jsol]*d0p[jsol][jsol]/(d0[jsol]*d0[jsol]))))*gradN[i]*mt.m_k[jsol]*z[jsol]*penalty + H[i]*H[j]*phif*vdzdc[isol][jsol]*mt.m_dkdc[jsol][jsol]*mt.m_c[jsol])*detJ;
                                }
                                else
                                {
                                    kcc += (((dDdc[ksol][jsol]*mt.m_k[ksol] + D[ksol]*mt.m_dkdc[ksol][jsol])*flux[ksol])*gradN[i]*H[j]*z[ksol]*penalty - (D[ksol]*fpt.m_w*H[j]*mt.m_c[ksol]*d0p[ksol][jsol]/(d0[ksol]*d0[ksol]))*gradN[i]*mt.m_k[ksol]*z[ksol]*penalty + H[i]*H[j]*phif*vdzdc[isol][jsol]*mt.m_dkdc[ksol][jsol]*mt.m_c[ksol])*detJ;
                                }
                            }
                        }
                        ke[i7+7+isol][j7+7+jsol] += kcc;
                    }
                    ke[i7  ][j7+7+isol] += kuc.x;
                    ke[i7+1][j7+7+isol] += kuc.y;
                    ke[i7+2][j7+7+isol] += kuc.z;
                    ke[i7+3][j7+7+isol] += kwc.x;
                    ke[i7+4][j7+7+isol] += kwc.y;
                    ke[i7+5][j7+7+isol] += kwc.z;
                    ke[i7+6][j7+7+isol] += kJc;
                    ke[i7+7+isol][j7  ] += kcu.x;
                    ke[i7+7+isol][j7+1] += kcu.y;
                    ke[i7+7+isol][j7+2] += kcu.z;
                    ke[i7+7+isol][j7+3] += kcw.x;
                    ke[i7+7+isol][j7+4] += kcw.y;
                    ke[i7+7+isol][j7+5] += kcw.z;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::StiffnessMatrix(FELinearSystem& LS)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
    int nsol = m_pMat->Solutes();
    int ndpn = 7 + nsol;
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        if (el.isActive()) {
            // element stiffness matrix
            FEElementMatrix ke(el);
            
            // create the element's stiffness matrix
            int ndof = ndpn*el.Nodes();
            ke.resize(ndof, ndof);
            ke.zero();
            
            // calculate material stiffness
            ElementStiffness(el, ke);
            
            // get the element's LM vector
            vector<int> lm;
            UnpackLM(el, lm);
            ke.SetIndices(lm);
            
            // assemble element matrix in global stiffness matrix
            LS.Assemble(ke);
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::MassMatrix(FELinearSystem& LS)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
    const int nsol = m_pMat->Solutes();
    const int ndpn = 7 + nsol;
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        if (el.isActive()) {
            
            FEElementMatrix ke(el);
            
            // create the element's stiffness matrix
            int ndof = ndpn*el.Nodes();
            ke.resize(ndof, ndof);
            ke.zero();
            
            // calculate inertial stiffness
            ElementMassMatrix(el, ke);
            
            // get the element's LM vector
            vector<int> lm;
            UnpackLM(el, lm);
            ke.SetIndices(lm);
            
            // assemble element matrix in global stiffness matrix
            LS.Assemble(ke);
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
    FEBiphasicFSI* pme = dynamic_cast<FEBiphasicFSI*>(GetMaterial()); assert(pme);
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
    const int nsol = m_pMat->Solutes();
    const int ndpn = 7 + nsol;
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        if (el.isActive()) {
            
            // element stiffness matrix
            FEElementMatrix ke(el);
            
            // create the element's stiffness matrix
            int ndof = ndpn*el.Nodes();
            ke.resize(ndof, ndof);
            ke.zero();
            
            // calculate inertial stiffness
            ElementBodyForceStiffness(bf, el, ke);
            
            // get the element's LM vector
            vector<int> lm;
            UnpackLM(el, lm);
            ke.SetIndices(lm);
            
            // assemble element matrix in global stiffness matrix
            LS.Assemble(ke);
        }
    }
}

//-----------------------------------------------------------------------------
//! calculates element inertial stiffness matrix
void FEMultiphasicFSIDomain3D::ElementMassMatrix(FESolidElement& el, matrix& ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i7, j, j7, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    double dtrans = m_btrans ? 1 : m_sseps;
    
    const int nsol = m_pMat->Solutes();
    const int ndpn = 7 + nsol;
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    vector<mat3d> gradgradN(neln);
    
    double *H;
    double *Gr, *Gs, *Gt;
    vec3d g[3], dg[3][3];
    
    // jacobian
    double Ji[3][3], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    double dt = tp.timeIncrement;
    double a = tp.gamma/(tp.beta*dt);
    double b = tp.alpham/(tp.alphaf*tp.beta*dt*dt);
    double c = tp.alpham/(tp.alphaf*tp.gamma*dt);
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjact(el, Ji, n, tp.alphaf)*gw[n]*tp.alphaf;
        
        ContraBaseVectors(el, n, g, tp.alphaf);
        ContraBaseVectorDerivatives(el, n, dg, tp.alphaf);
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // get the shape function derivatives
        double* Grr = el.Grr(n); double* Grs = el.Grs(n); double* Grt = el.Grt(n);
        double* Gsr = el.Gsr(n); double* Gss = el.Gss(n); double* Gst = el.Gst(n);
        double* Gtr = el.Gtr(n); double* Gts = el.Gts(n); double* Gtt = el.Gtt(n);
        
        // setup the material point
        // NOTE: deformation gradient and determinant have already been evaluated in the stress routine
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& et = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFSIMaterialPoint& fpt = *(mp.ExtractData<FEFSIMaterialPoint>());
        FEBiphasicFSIMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
        double Jf = 1 + pt.m_ef;
        
        double denss = m_pMat->SolidDensity(mp);
        double densTf = m_pMat->TrueFluidDensity(mp);
        
        // Jsdot/Js = div(vs)
        double dJsoJ = fpt.m_Jdot/et.m_J;
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate spatial gradgrad of shape functions
        for (i=0; i<neln; ++i) {
            gradgradN[i] = (((dg[0][0] & g[0]) + (dg[0][1] & g[1]) + (dg[0][2] & g[2]))*Gr[i]
                            + ((dg[1][0] & g[0]) + (dg[1][1] & g[1]) + (dg[1][2] & g[2]))*Gs[i]
                            + ((dg[2][0] & g[0]) + (dg[2][1] & g[1]) + (dg[2][2] & g[2]))*Gt[i]
                            + (g[0] & g[0])*Grr[i] + (g[0] & g[1])*Gsr[i] + (g[0] & g[2])*Gtr[i]
                            + (g[1] & g[0])*Grs[i] + (g[1] & g[1])*Gss[i] + (g[1] & g[2] )*Gts[i]
                            + (g[2] & g[0])*Grt[i] + (g[2] & g[1])*Gst[i] + (g[2] & g[2])*Gtt[i]);
        }
        
        double phif = m_pMat->Porosity(mp);
        double phis = m_pMat->SolidVolumeFrac(mp);
        vec3d gradphif = m_pMat->gradPorosity(mp);
        
        // evaluate stiffness matrix
        for (i=0, i7=0; i<neln; ++i, i7 += ndpn)
        {
            for (j=0, j7 = 0; j<neln; ++j, j7 += ndpn)
            {
                
                mat3d Kuu = ((mat3dd(b*H[j]*dtrans))*denss*H[i] + (((et.m_a*phis)&gradN[j]) + mat3dd(b*phif*H[j]*dtrans) - ((fpt.m_w&gradN[j]) * (-1.0/phif*dJsoJ + a*dtrans) - (fpt.m_w&(et.m_L.transpose()*gradN[j]*dtrans)))*phis/phif + mat3dd(gradN[j]*fpt.m_w*a*dtrans) - ((bpt.m_Lw*fpt.m_w)&gradN[j])*phis/(phif*phif) + ((fpt.m_w&gradN[j])*((gradphif*(phis+1.0)/phif)*fpt.m_w) - (fpt.m_w&(gradgradN[j].transpose()*fpt.m_w)))*phis/(phif*phif) - pt.m_Lf*(mat3dd(gradN[j]*fpt.m_w)) - (pt.m_aft&gradN[j])*phis)*(-H[i]*densTf*phis/phif))*detJ; //mixture 2
                
                mat3d Kuw = ((mat3dd(c*dtrans-phis/phif*dJsoJ-(gradphif*fpt.m_w)/(phif*phif))+pt.m_Lf)*H[j] + mat3dd(gradN[j]*fpt.m_w)/phif)*(-H[i]*densTf/phif*phis*detJ); //mixture 2
                
                vec3d kuJ = pt.m_aft*(densTf/Jf*H[i]*H[j]*phis*detJ); //mixture 2
                mat3d Kwu = (((et.m_a&gradN[j])*phis + mat3dd(b*phif*H[j]*dtrans) - ((fpt.m_w&gradN[j]) * (-1.0/phif*dJsoJ + a*dtrans) - (fpt.m_w&(et.m_L.transpose()*gradN[j]*dtrans)))*phis/phif + mat3dd(gradN[j]*fpt.m_w*a*dtrans) - ((bpt.m_Lw*fpt.m_w)&gradN[j])*phis/(phif*phif) + ((fpt.m_w&gradN[j])*((gradphif*(phis+1.0)/phif)*fpt.m_w) - (fpt.m_w&(gradgradN[j].transpose()*fpt.m_w)))*phis/(phif*phif) - pt.m_Lf*mat3dd(gradN[j]*fpt.m_w))*(H[i]*densTf/phif) + (pt.m_aft&gradN[j])*H[i]*densTf*(1.0-phis/phif))*detJ;
                mat3d Kww = ((mat3dd(c*dtrans-phis/phif*dJsoJ-(gradphif*fpt.m_w)/(phif*phif))+pt.m_Lf)*H[j] + mat3dd(gradN[j]*fpt.m_w)/phif)*(H[i]*densTf/phif*detJ);
                vec3d kwJ = pt.m_aft*(-densTf/Jf*H[i]*H[j]*detJ);
                
                ke[i7+0][j7  ] += Kuu(0,0); ke[i7+0][j7+1] += Kuu(0,1); ke[i7+0][j7+2] += Kuu(0,2);
                ke[i7+1][j7  ] += Kuu(1,0); ke[i7+1][j7+1] += Kuu(1,1); ke[i7+1][j7+2] += Kuu(1,2);
                ke[i7+2][j7  ] += Kuu(2,0); ke[i7+2][j7+1] += Kuu(2,1); ke[i7+2][j7+2] += Kuu(2,2);
                
                ke[i7+0][j7+3] += Kuw(0,0); ke[i7+0][j7+4] += Kuw(0,1); ke[i7+0][j7+5] += Kuw(0,2);
                ke[i7+1][j7+3] += Kuw(1,0); ke[i7+1][j7+4] += Kuw(1,1); ke[i7+1][j7+5] += Kuw(1,2);
                ke[i7+2][j7+3] += Kuw(2,0); ke[i7+2][j7+4] += Kuw(2,1); ke[i7+2][j7+5] += Kuw(2,2);
                
                ke[i7+3][j7  ] += Kwu(0,0); ke[i7+3][j7+1] += Kwu(0,1); ke[i7+3][j7+2] += Kwu(0,2);
                ke[i7+4][j7  ] += Kwu(1,0); ke[i7+4][j7+1] += Kwu(1,1); ke[i7+4][j7+2] += Kwu(1,2);
                ke[i7+5][j7  ] += Kwu(2,0); ke[i7+5][j7+1] += Kwu(2,1); ke[i7+5][j7+2] += Kwu(2,2);
                
                ke[i7+3][j7+3] += Kww(0,0); ke[i7+3][j7+4] += Kww(0,1); ke[i7+3][j7+5] += Kww(0,2);
                ke[i7+4][j7+3] += Kww(1,0); ke[i7+4][j7+4] += Kww(1,1); ke[i7+4][j7+5] += Kww(1,2);
                ke[i7+5][j7+3] += Kww(2,0); ke[i7+5][j7+4] += Kww(2,1); ke[i7+5][j7+5] += Kww(2,2);
                
                ke[i7+0][j7+6] += kuJ.x;
                ke[i7+1][j7+6] += kuJ.y;
                ke[i7+2][j7+6] += kuJ.z;
                ke[i7+3][j7+6] += kwJ.x;
                ke[i7+4][j7+6] += kwJ.y;
                ke[i7+5][j7+6] += kwJ.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::Update(const FETimeInfo& tp)
{
    bool berr = false;
    int NE = (int) m_Elem.size();
#pragma omp parallel for shared(NE, berr)
    for (int i=0; i<NE; ++i)
    {
        try
        {
            FESolidElement& el = Element(i);
            if (el.isActive())
            {
                UpdateElementStress(i, tp);
            }
        }
        catch (NegativeJacobian e)
        {
#pragma omp critical
            {
                // reset the logfile mode
                berr = true;
                if (e.DoOutput()) feLogError(e.what());
            }
        }
    }
    
    if (berr) throw NegativeJacobianDetected();
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
void FEMultiphasicFSIDomain3D::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double alphaf = tp.alphaf;
    double alpham = tp.alpham;
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    double dtrans = m_btrans ? 1 : m_sseps;
    
    // get the solid element
    FESolidElement& el = m_Elem[iel];
    
    // get the number of integration points
    int nint = el.GaussPoints();
    
    // number of nodes
    int neln = el.Nodes();
    
    // number of solutes
    const int nsol = m_pMat->Solutes();
    
    // nodal coordinates
    const int NELN = FEElement::MAX_NODES;
    vec3d r0[NELN], r[NELN];
    vec3d vs[NELN];
    vec3d a[NELN];
    vec3d w[NELN];
    vec3d aw[NELN];
    double e[NELN];
    double ae[NELN];
    vector< vector<double> > ct(nsol, vector<double>(NELN));
    vector< vector<double> > act(nsol, vector<double>(NELN));
    vector<int> sid(nsol);
    for (int j=0; j<nsol; ++j) sid[j] = m_pMat->GetSolute(j)->GetSoluteDOF();
    for (int j=0; j<neln; ++j) {
        FENode& node = m_pMesh->Node(el.m_node[j]);
        r0[j] = node.m_r0;
        r[j]  = node.m_rt*alphaf + node.m_rp*(1-alphaf);
        vs[j] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2])*alphaf + node.m_vp*(1-alphaf);
        w[j]  = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*alphaf + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1-alphaf);
        e[j]  = node.get(m_dofEF)*alphaf + node.get_prev(m_dofEF)*(1-alphaf);
        a[j]  = node.m_at*alpham + node.m_ap*(1-alpham);
        aw[j] = node.get_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2])*alpham + node.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2])*(1-alpham);
        ae[j] = node.get(m_dofAEF)*alpham + node.get_prev(m_dofAEF)*(1-alpham);
        for (int k=0; k<nsol; ++k) {
            ct[k][j] = node.get(m_dofC + sid[k])*alphaf + node.get_prev(m_dofC + sid[k])*(1-alphaf);
            act[k][j] = node.get(m_dofAC + sid[k])*alpham + node.get_prev(m_dofAC + sid[k])*(1-alpham);
        }
    }
    
    // loop over the integration points and update
    // velocity, velocity gradient, acceleration
    // stress and pressure at the integration point
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEFSIMaterialPoint& ft = *(mp.ExtractData<FEFSIMaterialPoint>());
        FEBiphasicFSIMaterialPoint& bt = *(mp.ExtractData<FEBiphasicFSIMaterialPoint>());
        FEMultiphasicFSIMaterialPoint& spt = *(mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
        
        // elastic material point data
        mp.m_r0 = el.Evaluate(r0, n);
        mp.m_rt = el.Evaluate(r, n);
        mat3d Ft, Fp;
        double Jt, Jp;
        Jt = defgrad(el, Ft, n);
        Jp = defgradp(el, Fp, n);
        ept.m_F = Ft*alphaf + Fp*(1 - alphaf);
        ept.m_J = ept.m_F.det();
        mat3d Fi = ept.m_F.inverse();
        ept.m_L = (Ft - Fp)*Fi*(dtrans / dt);
        ept.m_v = m_btrans ? el.Evaluate(vs, n) : vec3d(0, 0, 0);
        ept.m_a = m_btrans ? el.Evaluate(a, n) : vec3d(0, 0, 0);
        //NOTE: Made these new function. Converts from GradJ to gradJ
        vec3d GradJ = FESolidDomain::GradJ(el,n);
        vec3d GradJp = FESolidDomain::GradJp(el, n);
        
        //calculate gradJ gradphif
        bt.m_gradJ = Fi.transpose()*(GradJ*alphaf + GradJp*(1-alphaf));
        ept.m_gradJ = bt.m_gradJ;
        
        // FSI material point data
        ft.m_w = el.Evaluate(w, n);
        ft.m_Jdot = (Jt - Jp) / dt*dtrans;
        ft.m_aw = el.Evaluate(aw, n)*dtrans;
        
        for (int isol=0; isol < nsol; ++isol) {
            spt.m_c[isol] = el.Evaluate(ct[isol], n);
            vec3d Gradc = Gradient(el, ct[isol], n);
            spt.m_gradc[isol] = Fi.transpose()*Gradc;
            spt.m_cdot[isol] = el.Evaluate(act[isol], n)*dtrans;
        }
        
        bt.m_phi0 = m_pMat->SolidReferentialVolumeFraction(mp);
        m_pMat->PartitionCoefficientFunctions(mp, spt.m_k, spt.m_dkdJ, spt.m_dkdc);
        
        // fluid material point data
        double phif = m_pMat->Porosity(mp);
        double phis = m_pMat->SolidVolumeFrac(mp);
        vec3d gradphif = m_pMat->gradPorosity(mp);
        pt.m_efdot = el.Evaluate(ae, n)*dtrans;
        pt.m_vft = ept.m_v + ft.m_w/phif;
        mat3d Gradw = Gradient(el, w, n);
        bt.m_Lw = Gradw*Fi;
        pt.m_Lf = ept.m_L + bt.m_Lw/phif - (ft.m_w & gradphif)/(phif*phif);
        pt.m_ef = el.Evaluate(e, n);
        vec3d Gradef = Gradient(el, e, n);
        pt.m_gradef = Fi.transpose()*Gradef;
        
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        
        // fluid acceleration
        pt.m_aft = (ft.m_aw + ept.m_a*phif - ft.m_w*ft.m_Jdot*phis/phif/ept.m_J + pt.m_Lf*ft.m_w)/phif;
        
        spt.m_pe = m_pMat->Fluid()->Pressure(mp);
        
        // calculate the solute flux and actual concentration
        for (int isol=0; isol < nsol; ++isol)
        {
            spt.m_j[isol] = m_pMat->SoluteFlux(mp, isol);
            spt.m_ca[isol] = m_pMat->ConcentrationActual(mp, isol);
        }
        
        // calculate the fluid pressure
        pt.m_pf = m_pMat->PressureActual(mp);
        
        spt.m_psi = m_pMat->ElectricPotential(mp);
        spt.m_Ie = m_pMat->CurrentDensity(mp);
        spt.m_cF = m_pMat->FixedChargeDensity(mp);
        
        // calculate the solid stress at this material point
        pt.m_sf = m_pMat->Fluid()->GetViscous()->Stress(mp); //Old Nat BC
        ft.m_ss = m_pMat->Solid()->Stress(mp) - m_pMat->Fluid()->GetViscous()->Stress(mp)*phis; //Old NatBC
         //Old Nat BC
        for (int isol=0; isol<nsol; ++isol)
            ft.m_ss += -mat3dd(1.0)*R*T*osmc*spt.m_ca[isol];
        
        
        //pt.m_sf = m_pMat->Fluid()->GetViscous()->Stress(mp); //New Nat BC
        //ft.m_ss = -mat3dd(1.0)*pt.m_pf + m_pMat->Solid()->Stress(mp) - m_pMat->Fluid()->GetViscous()->Stress(mp)*phis; //New NatBC
        
        // calculate the mixture stress at this material point
        ept.m_s =  -mat3dd(1.0)*spt.m_pe + ft.m_ss + pt.m_sf;
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::InertialForces(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
    
    const int nsol = m_pMat->Solutes();
    const int ndpn = 7+nsol;
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // get the element
        FESolidElement& el = m_Elem[i];
        
        if (el.isActive()) {
            // element force vector
            vector<double> fe;
            vector<int> lm;
            
            // get the element force vector and initialize it to zero
            int ndof = ndpn*el.Nodes();
            fe.assign(ndof, 0);
            
            // calculate internal force vector
            ElementInertialForce(el, fe);
            
            // get the element's LM vector
            UnpackLM(el, lm);
            
            // assemble element 'fe'-vector into global R vector
            R.Assemble(el.m_node, lm, fe);
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::ElementInertialForce(FESolidElement& el, vector<double>& fe)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, n;
    
    // jacobian determinant
    double detJ;
    
    const double* H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    const int ndpn = 7 + nsol;
    
    double*    gw = el.GaussWeights();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
        double densTf = m_pMat->TrueFluidDensity(mp);
        double densTs = m_pMat->TrueSolidDensity(mp);
        double phis = m_pMat->SolidVolumeFrac(mp);
        
        // calculate the jacobian
        detJ = detJt(el, n, tp.alphaf)*gw[n];
        
        H = el.H(n);
        
        for (i=0; i<neln; ++i)
        {
            vec3d f = pt.m_aft*(densTf*H[i]*detJ);
            
            vec3d fs = (-pt.m_aft*densTf + ept.m_a*densTs)*H[i]*phis*detJ; //mixture 2
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[ndpn*i+0] -= fs.x;
            fe[ndpn*i+1] -= fs.y;
            fe[ndpn*i+2] -= fs.z;
            fe[ndpn*i+3] -= f.x;
            fe[ndpn*i+4] -= f.y;
            fe[ndpn*i+5] -= f.z;
        }
    }
}

//-----------------------------------------------------------------------------
void FEMultiphasicFSIDomain3D::Serialize(DumpStream& ar)
{
    FESolidDomain::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_sseps;
    ar & m_pMat;
    ar & m_dofU & m_dofV & m_dofW & m_dofAW;
    ar & m_dofSU & m_dofR;
    ar & m_dof;
    ar & m_dofEF & m_dofAEF & m_dofC & m_dofAC;
}
