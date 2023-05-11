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
#include "FEFluidSolutesDomain3D.h"
#include "FEFluidSolver.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioFluidSolutes.h"
#include <FECore/FELinearSystem.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidSolutesDomain3D, FESolidDomain)
    ADD_PARAMETER(m_newFormulation   , "new_formula");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEFluidSolutesDomain3D::FEFluidSolutesDomain3D(FEModel* pfem) : FESolidDomain(pfem), FEFluidDomain(pfem), m_dofW(pfem), m_dofAW(pfem), m_dof(pfem)
{
    m_pMat = 0;
    m_btrans = true;
    m_newFormulation = true;
    
    // TODO: Can this be done in Init, since  there is no error checking
    if (pfem)
    {
        m_dofW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY));
        m_dofAW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_ACCELERATION));

        m_dofEF = pfem->GetDOFIndex("ef");
        m_dofAEF = pfem->GetDOFIndex("aef");
        m_dofC = pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
        m_dofAC = pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION_TDERIV), 0);

        // list the degrees of freedom
        // (This allows the FEDomain base class to handle several tasks such as UnpackLM)
        m_dof.AddDof(m_dofW[0]);
        m_dof.AddDof(m_dofW[1]);
        m_dof.AddDof(m_dofW[2]);
        m_dof.AddDof(m_dofEF);
    }
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEFluidSolutesDomain3D& FEFluidSolutesDomain3D::operator = (FEFluidSolutesDomain3D& d)
{
    m_Elem = d.m_Elem;
    m_pMesh = d.m_pMesh;
    return (*this);
}

//-----------------------------------------------------------------------------
//! get the total dofs
const FEDofList& FEFluidSolutesDomain3D::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
//! Assign material
void FEFluidSolutesDomain3D::SetMaterial(FEMaterial* pmat)
{
    FEDomain::SetMaterial(pmat);
    if (pmat)
    {
        m_pMat = dynamic_cast<FEFluidSolutes*>(pmat);
        assert(m_pMat);
    }
    else m_pMat = 0;
}

//-----------------------------------------------------------------------------
bool FEFluidSolutesDomain3D::Init()
{
    // initialize base class
    if (FESolidDomain::Init() == false) return false;
    
    const int nsol = m_pMat->Solutes();
    
    // set the active degrees of freedom list
    FEDofList dofs = GetDOFList();
    for (int i=0; i<nsol; ++i)
    {
        int m = m_pMat->GetSolute(i)->GetSoluteDOF();
        dofs.AddDof(m_dofC + m);
    }
    m_dof = dofs;
    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::Serialize(DumpStream& ar)
{
    FESolidDomain::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_pMat;
    ar & m_dofW & m_dofAW & m_dof;
    ar & m_dofEF & m_dofAEF & m_dofC & m_dofAC;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::Reset()
{
    // reset base class
    FESolidDomain::Reset();
    
    const int nsol = m_pMat->Solutes();
    
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
            FEFluidSolutesMaterialPoint& ps = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
            
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
            
            for (int j=0; j<m_pMat->Reactions(); ++j)
                m_pMat->GetReaction(j)->ResetElementData(mp);
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::Activate()
{
    const int nsol = m_pMat->Solutes();
    
    for (int i=0; i<Nodes(); ++i)
    {
        FENode& node = Node(i);
        if (node.HasFlags(FENode::EXCLUDE) == false)
        {
/*          if (node.m_rid < 0)
            {
                node.set_active(m_dofU[0]);
                node.set_active(m_dofU[1]);
                node.set_active(m_dofU[2]);
            }
*/            node.set_active(m_dofW[0]);
            node.set_active(m_dofW[1]);
            node.set_active(m_dofW[2]);
            node.set_active(m_dofEF);
            for (int isol=0; isol<nsol; ++isol)
                node.set_active(m_dofC + m_pMat->GetSolute(isol)->GetSoluteDOF());
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::InitMaterialPoints()
{
    const int nsol = m_pMat->Solutes();
    FEMesh& m = *GetMesh();
    
    const int NE = FEElement::MAX_NODES;
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
            for (int isol = 0; isol<nsol; ++isol)
                c0[isol][i] = ni.get(m_dofC + sid[isol]);
        }

        // get the number of integration points
        int nint = el.GaussPoints();
        
        // loop over the integration points
        for (int n = 0; n<nint; ++n)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(n);
            FEFluidSolutesMaterialPoint& ps = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
            
            // initialize solutes
            ps.m_nsol = nsol;
            
            // initialize effective solute concentrations
            for (int isol = 0; isol<nsol; ++isol) {
                ps.m_c[isol] = el.Evaluate(c0[isol], n);
                ps.m_ca[isol] = m_pMat->ConcentrationActual(mp, isol);
                ps.m_gradc[isol] = gradient(el, c0[isol], n);
            }
            
            ps.m_psi = m_pMat->ElectricPotential(mp);
            ps.m_Ie = m_pMat->CurrentDensity(mp);
            
            for (int isol = 0; isol<nsol; ++isol)
                ps.m_j[isol] = m_newFormulation ? m_pMat->SoluteDiffusiveFlux(mp, isol) : m_pMat->SoluteFlux(mp, isol);
        }
    }
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEFluidSolutesDomain3D::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    const int NE = FEElement::MAX_NODES;
    vec3d x0[NE], r0, v;
    FEMesh& m = *GetMesh();
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FESolidElement& el = m_Elem[i];
        int neln = el.Nodes();
        for (int i=0; i<neln; ++i)
        {
            x0[i] = m.Node(el.m_node[i]).m_r0;
        }
        
        int n = el.GaussPoints();
        for (int j=0; j<n; ++j)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(j);
            FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
            pt.m_r0 = el.Evaluate(x0, j);
            mp.m_rt = mp.m_r0;
            
            if (pt.m_ef <= -1) {
                throw NegativeJacobianDetected();
            }
            
            // reset chemical reaction element data
            for (int j=0; j<m_pMat->Reactions(); ++j)
                m_pMat->GetReaction(j)->InitializeElementData(mp);
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::InternalForces(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
    
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
#pragma omp parallel for shared (NE)
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
        if (m_newFormulation) ElementInternalForceNew(el, fe);
        else ElementInternalForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEFluidSolutesDomain3D::ElementInternalForceNew(FESolidElement& el, vector<double>& fe)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJ;
    
    mat3ds sv;
    vec3d gradep;
    
    const double *H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    const int nreact = m_pMat->Reactions();

    double dt = tp.timeIncrement;
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double*    gw = el.GaussWeights();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());

        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        // get the viscous stress tensor for this integration point
        sv = m_pMat->Fluid()->GetViscous()->Stress(mp);
        // get the gradient of the elastic pressure
        gradep = pt.m_gradef*m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        
        // Miscellaneous constants
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double penalty = m_pMat->m_penalty;
        
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
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
        {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        // Jdot/J
        double dJoJ = pt.m_efdot/(pt.m_ef+1);
        
        double dms = m_pMat->m_diffMtmSupp;
        
        vector<int> z(nsol);
        vector<vec3d> jd(nsol);
        vec3d je(0,0,0);
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            jd[isol] = m_pMat->SoluteDiffusiveFlux(mp,isol);
            je += jd[isol]*z[isol];
        }
        
        vector<double> dkdt(nsol,0);
        vector<vec3d> gradk(nsol,vec3d(0,0,0));
        for (int isol=0; isol<nsol; ++isol)
        {
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                dkdt[isol] += spt.m_dkdc[isol][jsol]*spt.m_cdot[jsol];
                gradk[isol] += spt.m_gradc[jsol]*spt.m_dkdc[isol][jsol];
            }
        }
        
        for (i=0; i<neln; ++i)
        {
            vec3d fs = sv*gradN[i] + gradep*H[i];
            for (int isol=0; isol<nsol; ++isol)
                fs += spt.m_gradc[isol]*(R*T*spt.m_k[isol]*H[i]*dms); //fluid mtm bal only
            double fJ = (dJoJ+phiwhat)*H[i] + gradN[i]*pt.m_vft;
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[ndpn*i  ] -= fs.x*detJ;
            fe[ndpn*i+1] -= fs.y*detJ;
            fe[ndpn*i+2] -= fs.z*detJ;
            fe[ndpn*i+3] -= fJ*detJ;
            for (int isol=0; isol<nsol; ++isol)
            {
                double fc = (jd[isol]+je*penalty)*gradN[i] - H[i]*((spt.m_cdot[isol]+spt.m_gradc[isol]*pt.m_vft)*spt.m_k[isol] + spt.m_c[isol]*(dkdt[isol]+gradk[isol]*pt.m_vft) - chat[isol] + spt.m_ca[isol]/(1+pt.m_ef)*(pt.m_efdot + pt.m_gradef*pt.m_vft));
                fe[ndpn*i+4+isol] -= fc*detJ;
            }
        }
    }
}

void FEFluidSolutesDomain3D::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJ;
    
    mat3ds sv;
    vec3d gradep;
    
    const double *H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    const int nreact = m_pMat->Reactions();

    double dt = tp.timeIncrement;
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double*    gw = el.GaussWeights();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());

        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        // get the viscous stress tensor for this integration point
        sv = m_pMat->Fluid()->GetViscous()->Stress(mp);
        // get the gradient of the elastic pressure
        gradep = pt.m_gradef*m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        
        // Miscellaneous constants
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double penalty = m_pMat->m_penalty;
        
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
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
        {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        // Jdot/J
        double dJoJ = pt.m_efdot/(pt.m_ef+1);
        
        double dms = m_pMat->m_diffMtmSupp;
        
        vector<int> z(nsol);
        vec3d je(0,0,0);
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            je += spt.m_j[isol]*z[isol];
        }
        
        vector<double> dkdt(nsol,0);
        for (int isol=0; isol<nsol; ++isol)
        {
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                dkdt[isol] += spt.m_dkdc[isol][jsol]*spt.m_cdot[jsol];
            }
        }
        
        for (i=0; i<neln; ++i)
        {
            vec3d fs = sv*gradN[i] + gradep*H[i];
            for (int isol=0; isol<nsol; ++isol)
                fs += spt.m_gradc[isol]*(R*T*spt.m_k[isol]*H[i]*dms); //fluid mtm bal only
            double fJ = (dJoJ+phiwhat)*H[i] + gradN[i]*pt.m_vft;
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[ndpn*i  ] -= fs.x*detJ;
            fe[ndpn*i+1] -= fs.y*detJ;
            fe[ndpn*i+2] -= fs.z*detJ;
            fe[ndpn*i+3] -= fJ*detJ;
            for (int isol=0; isol<nsol; ++isol)
            {
                double fc = (spt.m_j[isol]+je*penalty)*gradN[i] - H[i]*(spt.m_cdot[isol]*spt.m_k[isol] + spt.m_c[isol]*dkdt[isol] - chat[isol]);
                fe[ndpn*i+4+isol] -= fc*detJ;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
    
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = ndpn*el.Nodes();
        fe.assign(ndof, 0);
        
        // apply body forces
        if (m_newFormulation) ElementBodyForceNew(BF, el, fe);
        else ElementBodyForce(BF, el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEFluidSolutesDomain3D::ElementBodyForceNew(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
{
    double R = m_pMat->m_Rgas;
    double T = m_pMat->m_Tabs;
    
    // jacobian
    double Ji[3][3], detJ;
    const double *H, *Gr, *Gs, *Gt;
    double* gw = el.GaussWeights();
    vec3d f;
    
    // number of nodes
    int neln = el.Nodes();
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
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
        FEFluidSolutesMaterialPoint& spt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
        double dens = m_pMat->Fluid()->Density(mp);
        
        pt.m_r0 = el.Evaluate(r0, n);
        
        detJ = invjac0(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        // get the force
        f = BF.force(mp);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        double dms = m_pMat->m_diffMtmSupp;
        double penalty = m_pMat->m_penalty;
        vector<vec3d> jb(nsol,vec3d(0,0,0));
        vec3d je(0,0,0);
        
        // evaluate sedimentation fluxes and effective flux contribution
        for (int isol=0; isol<nsol; ++isol) {
            double M = m_pMat->GetSolute(isol)->MolarMass();
            // get the sedimentation coefficient
            double s = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp)*M/(R*T);
            // get the sedimentation flux
            jb[isol] = f*(s*spt.m_ca[isol]);
            // get the charge number
            double z = m_pMat->GetSolute(isol)->ChargeNumber();
            je += jb[isol]*z;
            dens += M*spt.m_ca[isol]*dms;
        }
        
        // evaluate spatial gradient of shape functions
        for (int i=0; i<neln; ++i)
        {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        for (int i=0; i<neln; ++i)
        {
            vec3d fs = f*H[i]*dens;
            
            fe[ndpn*i  ] -= fs.x*detJ;
            fe[ndpn*i+1] -= fs.y*detJ;
            fe[ndpn*i+2] -= fs.z*detJ;
            for (int isol=0; isol<nsol; ++isol)
            {
                double fc = -gradN[i]*(jb[isol] + je*penalty);
                fe[ndpn*i+4+isol] -= fc*detJ;
            }
        }
    }
}

void FEFluidSolutesDomain3D::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
{
    double R = m_pMat->m_Rgas;
    double T = m_pMat->m_Tabs;
    
    // jacobian
    double Ji[3][3], detJ;
    const double *H, *Gr, *Gs, *Gt;
    double* gw = el.GaussWeights();
    vec3d f;
    
    // number of nodes
    int neln = el.Nodes();
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
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
        FEFluidSolutesMaterialPoint& spt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
        double dens = m_pMat->Fluid()->Density(mp);
        
        pt.m_r0 = el.Evaluate(r0, n);
        
        detJ = invjac0(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        // get the force
        f = BF.force(mp);
        
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
        
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            M[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
        }
        
        // evaluate spatial gradient of shape functions
        for (int i=0; i<neln; ++i)
        {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        for (int i=0; i<neln; ++i)
        {
            vec3d fs = f*H[i]*dens;
            for(int isol = 0; isol<nsol; ++isol)
            {
                fs += f*H[i]*spt.m_k[isol]*spt.m_c[isol]*M[isol]*dms;
            }
            
            fe[ndpn*i  ] -= fs.x*detJ;
            fe[ndpn*i+1] -= fs.y*detJ;
            fe[ndpn*i+2] -= fs.z*detJ;
            for (int isol=0; isol<nsol; ++isol)
            {
                double fc = -gradN[i]*f*d0[isol]*spt.m_k[isol]*M[isol]*spt.m_c[isol]/(R*T);
                for(int jsol = 0; jsol<nsol; ++jsol)
                {
                    fc += -gradN[i]*f*z[jsol]*penalty*d0[jsol]*spt.m_k[jsol]*M[jsol]*spt.m_c[jsol]/(R*T);
                }
                fe[ndpn*i+4+isol] -= fc*detJ;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEFluidSolutesDomain3D::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    // jacobian
    double Ji[3][3], detJ;
    const double *H, *Gr, *Gs, *Gt;
    double* gw = el.GaussWeights();
    
    // number of nodes
    int neln = el.Nodes();
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    vec3d f, k;
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
        
        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n]*tp.alphaf;
        
        H = el.H(n);
        
        double dens = m_pMat->Fluid()->Density(mp);
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double dms = m_pMat->m_diffMtmSupp;
        double penalty = m_pMat->m_penalty;
        vector<double> M(nsol);
        vector<int> z(nsol);
        vector<double> d0(nsol);
        vector<vector<double>> d0p(nsol, vector<double>(nsol));
        
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            M[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                d0p[isol][jsol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, jsol);
            }
        }
        
        // get the force
        f = BF.force(mp);
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (int i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        for (int i=0, i4=0; i<neln; ++i, i4 += ndpn) {
            for (int j=0, j4 = 0; j<neln; ++j, j4 += ndpn)
            {
                k = f*(-H[i]*H[j]*dens/(pt.m_ef+1)*detJ);
                ke[i4  ][j4+3] += k.x;
                ke[i4+1][j4+3] += k.y;
                ke[i4+2][j4+3] += k.z;
                
                for (int isol = 0; isol<nsol; ++isol)
                {
                    vec3d kvc = vec3d(0.0);
                    for(int jsol=0; jsol<nsol; ++jsol)
                    {
                        double kcc = 0;
                        if (isol == jsol)
                        {
                            kvc += f*(spt.m_dkdc[isol][isol]*spt.m_c[isol]+spt.m_k[isol])*M[isol]*H[i]*H[j];
                            kcc = -gradN[i]*f*(spt.m_dkdc[isol][isol]*d0[isol]*spt.m_c[isol]+spt.m_k[isol]*d0p[isol][isol]*spt.m_c[isol]+spt.m_k[isol]*d0[isol])*M[isol]*H[j]/(R*T);
                            for (int ksol=0; ksol<nsol; ++ksol)
                            {
                                if(isol==ksol)
                                {
                                    kcc += -gradN[i]*f*(spt.m_dkdc[isol][isol]*d0[isol]*spt.m_c[isol]+spt.m_k[isol]*d0p[isol][isol]*spt.m_c[isol]+spt.m_k[isol]*d0[isol])*M[isol]*z[isol]*penalty*H[j]/(R*T);
                                }
                                else
                                {
                                    kcc += -gradN[i]*f*(spt.m_dkdc[ksol][isol]*d0[ksol]+spt.m_k[ksol]*d0p[ksol][isol])*spt.m_c[ksol]*M[ksol]*z[ksol]*penalty*H[j]/(R*T);
                                }
                            }
                        }
                        else
                        {
                            kvc += f*spt.m_c[jsol]*spt.m_dkdc[jsol][isol]*M[jsol]*H[i]*H[j];
                            kcc = -gradN[i]*f*H[j]*(spt.m_dkdc[isol][jsol]*d0[isol]+spt.m_k[isol]*d0p[isol][jsol])*M[isol]*spt.m_c[isol]/(R*T);
                            for (int ksol=0; ksol<nsol; ++ksol)
                            {
                                if(jsol==ksol)
                                {
                                    kcc += -gradN[i]*f*(spt.m_dkdc[jsol][jsol]*d0[jsol]*spt.m_c[jsol]+spt.m_k[jsol]*d0p[jsol][jsol]*spt.m_c[jsol]+spt.m_k[jsol]*d0[jsol])*M[jsol]*z[jsol]*penalty*H[j]/(R*T);
                                }
                                else
                                {
                                    kcc += -gradN[i]*f*(spt.m_dkdc[ksol][jsol]*d0[ksol]+spt.m_k[ksol]*d0p[ksol][jsol])*spt.m_c[ksol]*M[ksol]*z[ksol]*penalty*H[j]/(R*T);
                                }
                            }
                        }
                        ke[i4+4+isol][j4+4+jsol] += kcc*detJ;
                    }
                    ke[i4  ][j4+4+isol] += kvc.x*detJ*dms;
                    ke[i4+1][j4+4+isol] += kvc.y*detJ*dms;
                    ke[i4+2][j4+4+isol] += kvc.z*detJ*dms;
                }
            }
        }
    }
}

void FEFluidSolutesDomain3D::ElementBodyForceStiffnessNew(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    // jacobian
    double Ji[3][3], detJ;
    const double *H, *Gr, *Gs, *Gt;
    double* gw = el.GaussWeights();
    
    // number of nodes
    int neln = el.Nodes();
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    vec3d f, k;
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
        
        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n]*tp.alphaf;
        
        H = el.H(n);
        
        double dens = m_pMat->Fluid()->Density(mp);
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double dms = m_pMat->m_diffMtmSupp;
        double penalty = m_pMat->m_penalty;
        vector<double> M(nsol);
        vector<int> z(nsol);
        vector<double> d0(nsol);
        vector<double> s(nsol);
        vector<vector<double>> d0p(nsol, vector<double>(nsol));
        
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            M[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            s[isol] = d0[isol]*M[isol]/(R*T);
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                d0p[isol][jsol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, jsol);
            }
        }
        
        // get the force
        f = BF.force(mp);
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (int i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        for (int i=0, i4=0; i<neln; ++i, i4 += ndpn) {
            for (int j=0, j4 = 0; j<neln; ++j, j4 += ndpn)
            {
                k = f*(-H[i]*H[j]*dens/(pt.m_ef+1)*detJ);
                ke[i4  ][j4+3] += k.x;
                ke[i4+1][j4+3] += k.y;
                ke[i4+2][j4+3] += k.z;
                
                for (int isol = 0; isol<nsol; ++isol)
                {
                    vec3d kvc = vec3d(0.0);
                    double kcc = (gradN[i]*f)*(s[isol]*H[j]);
                    for(int jsol=0; jsol<nsol; ++jsol)
                    {
                        if (isol == jsol)
                        {
                            kcc *= 2;
                            kvc += f*(spt.m_dkdc[isol][isol]*spt.m_c[isol]+spt.m_k[isol])*M[isol]*H[i]*H[j];
                        }
                        else
                        {
                            kvc += f*spt.m_c[jsol]*spt.m_dkdc[jsol][isol]*M[jsol]*H[i]*H[j];
                        }
                        ke[i4+4+isol][j4+4+jsol] += kcc*detJ;
                    }
                    ke[i4  ][j4+4+isol] += kvc.x*detJ*dms;
                    ke[i4+1][j4+4+isol] += kvc.y*detJ*dms;
                    ke[i4+2][j4+4+isol] += kvc.z*detJ*dms;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEFluidSolutesDomain3D::ElementStiffnessNew(FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i4, j, j4, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    const int ndpn = 4 + nsol;
    
    const int nreact = m_pMat->Reactions();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);
    
    double *H, *Gr, *Gs, *Gt;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjac0(el, Ji, n)*gw[n]*tp.alphaf;
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // setup the material point
        // NOTE: deformation gradient and determinant have already been evaluated in the stress routine
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
        double Jf = 1 + pt.m_ef;

        // get the tangents
        mat3ds svJ = m_pMat->Fluid()->GetViscous()->Tangent_Strain(mp);
        tens4ds cv = m_pMat->Fluid()->Tangent_RateOfDeformation(mp);
        double dep = m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        double d2ep = m_pMat->Fluid()->Tangent_Pressure_Strain_Strain(mp);
        // Jdot/J
        double dJoJ = pt.m_efdot/Jf;
        // Miscellaneous constants
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double dms = m_pMat->m_diffMtmSupp;
        double penalty = m_pMat->m_penalty;
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        double dodJ = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Strain(mp);
        vector<double> dodc(nsol);
        vector<double> M(nsol);
        vector<int> z(nsol);
        double kzc = 0;
        vector<double> d0(nsol);
        vector<vector<double>> d0p(nsol, vector<double>(nsol));
        
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            dodc[isol] = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp,isol);
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            M[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            kzc += z[isol]*spt.m_c[isol]*spt.m_k[isol];
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                d0p[isol][jsol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, jsol);
            }
        }
        
        //Get dk/dt (partial differential wrt time)
        vector<double> dkdt(nsol,0);
        vector<vec3d> gradk(nsol,vec3d(0,0,0));
        for (int isol=0; isol<nsol; ++isol)
        {
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                dkdt[isol] += spt.m_dkdc[isol][jsol]*spt.m_cdot[jsol];
                gradk[isol] += spt.m_gradc[jsol]*spt.m_dkdc[isol][jsol];
            }
        }
        
        // evaluate the chat
        vector<double> vbardzdc(nsol, 0.0);
        vector<vector<double>> vdzdc(nsol, vector<double>(nsol, 0.0));
        
        // chemical reactions
        for (i=0; i<nreact; ++i) {
            double vbar = m_pMat->GetReaction(i)->m_Vbar;
            for (int isol = 0; isol < nsol; ++isol)
            {
                vbardzdc[isol] += vbar*m_pMat->GetReaction(i)->Tangent_ReactionSupply_Concentration(mp,isol);
                for (int jsol = 0; jsol < nsol; ++jsol)
                {
                    double v = m_pMat->GetReaction(i)->m_v[isol];
                    vdzdc[isol][jsol] = v*m_pMat->GetReaction(i)->Tangent_ReactionSupply_Concentration(mp,jsol);
                }
            }
        }
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i4=0; i<neln; ++i, i4 += ndpn)
        {
            for (j=0, j4 = 0; j<neln; ++j, j4 += ndpn)
            {
                mat3d Kvv = vdotTdotv(gradN[i], cv, gradN[j]);
                vec3d kJv = (pt.m_gradef*(H[i]/Jf) + gradN[i])*H[j];
                vec3d kvJ = (svJ*gradN[i])*H[j] + (gradN[j]*dep+pt.m_gradef*H[j]*d2ep)*H[i];
                double kJJ = (H[j]*(ksi/dt - dJoJ) + gradN[j]*pt.m_vft)*H[i]/Jf;
                
                ke[i4  ][j4  ] += Kvv(0,0)*detJ;
                ke[i4  ][j4+1] += Kvv(0,1)*detJ;
                ke[i4  ][j4+2] += Kvv(0,2)*detJ;
                ke[i4  ][j4+3] += kvJ.x*detJ;
                
                ke[i4+1][j4  ] += Kvv(1,0)*detJ;
                ke[i4+1][j4+1] += Kvv(1,1)*detJ;
                ke[i4+1][j4+2] += Kvv(1,2)*detJ;
                ke[i4+1][j4+3] += kvJ.y*detJ;
                
                ke[i4+2][j4  ] += Kvv(2,0)*detJ;
                ke[i4+2][j4+1] += Kvv(2,1)*detJ;
                ke[i4+2][j4+2] += Kvv(2,2)*detJ;
                ke[i4+2][j4+3] += kvJ.z*detJ;
                
                ke[i4+3][j4  ] += kJv.x*detJ;
                ke[i4+3][j4+1] += kJv.y*detJ;
                ke[i4+3][j4+2] += kJv.z*detJ;
                ke[i4+3][j4+3] += kJJ*detJ;
                
                for (int isol=0; isol<nsol; ++isol) {
                    vec3d kcv = (pt.m_gradef*spt.m_ca[isol]/(1+pt.m_ef)+spt.m_gradc[isol]*spt.m_k[isol]
                                 +gradk[isol]*spt.m_c[isol])*(-H[i]*H[j]);
                    vec3d kvc(0,0,0);
                    double kJc = H[i]*H[j]*vbardzdc[isol]*spt.m_k[isol];
                    double kcJ = H[i]/(1+pt.m_ef)*(H[j]*(ksi/dt*spt.m_ca[isol] + spt.m_cdot[isol] + pt.m_vft*spt.m_gradc[isol])
                                                   +spt.m_ca[isol]*(pt.m_vft*gradN[j]));
                    int irow = i4+4+isol;
                    
                    for(int jsol=0; jsol<nsol; ++jsol)
                    {
                        kJc += H[i]*H[j]*vbardzdc[isol]*spt.m_dkdc[jsol][isol]*spt.m_c[jsol];
                        double kcc = H[i]*H[j]*vdzdc[isol][jsol] - z[jsol]*spt.m_k[jsol]*d0[jsol]*(gradN[i]*gradN[j]);
                        if (isol == jsol)
                        {
                            kvc += (spt.m_gradc[isol]*spt.m_dkdc[isol][isol]*H[j]+gradN[j]*spt.m_k[isol])*H[i]*R*T;
                            kcc -= H[i]*H[j]*(dJoJ + spt.m_dkdc[isol][isol]*(spt.m_cdot[isol]+ksi/dt*spt.m_c[isol]) + dkdt[isol] + ksi/dt*spt.m_k[isol] + (pt.m_vft*pt.m_gradef)) + gradN[j]*(pt.m_vft*H[i] + gradN[i]*(spt.m_k[isol]*d0[isol]));
                        }
                        else
                        {
                            kvc += spt.m_gradc[jsol]*spt.m_dkdc[jsol][isol]*H[i]*H[j]*R*T;
                        }
                        ke[irow][j4+4+jsol] += kcc*detJ;
                    }
                    ke[i4  ][j4+4+isol] += kvc.x*detJ*dms;
                    ke[i4+1][j4+4+isol] += kvc.y*detJ*dms;
                    ke[i4+2][j4+4+isol] += kvc.z*detJ*dms;
                    ke[i4+3][j4+4+isol] += kJc*detJ;
                    ke[irow][j4  ] += kcv.x*detJ;
                    ke[irow][j4+1] += kcv.y*detJ;
                    ke[irow][j4+2] += kcv.z*detJ;
                    ke[irow][j4+3] += kcJ*detJ;
                }
            }
        }
    }
}

void FEFluidSolutesDomain3D::ElementStiffness(FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i4, j, j4, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    const int ndpn = 4 + nsol;
    
    const int nreact = m_pMat->Reactions();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);
    
    double *H, *Gr, *Gs, *Gt;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjac0(el, Ji, n)*gw[n]*tp.alphaf;
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // setup the material point
        // NOTE: deformation gradient and determinant have already been evaluated in the stress routine
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
        double Jf = 1 + pt.m_ef;

        // get the tangents
        mat3ds svJ = m_pMat->Fluid()->GetViscous()->Tangent_Strain(mp);
        tens4ds cv = m_pMat->Fluid()->Tangent_RateOfDeformation(mp);
        double dep = m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        double d2ep = m_pMat->Fluid()->Tangent_Pressure_Strain_Strain(mp);
        // Jdot/J
        double dJoJ = pt.m_efdot/Jf;
        // Miscellaneous constants
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        double dms = m_pMat->m_diffMtmSupp;
        double penalty = m_pMat->m_penalty;
        double osmc = m_pMat->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
        double dodJ = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Strain(mp);
        vector<double> dodc(nsol);
        vector<double> M(nsol);
        vector<int> z(nsol);
        double kzc = 0;
        vector<double> d0(nsol);
        vector<vector<double>> d0p(nsol, vector<double>(nsol));
        
        for (int isol=0; isol<nsol; ++isol) {
            // get the charge number
            dodc[isol] = m_pMat->GetOsmoticCoefficient()->Tangent_OsmoticCoefficient_Concentration(mp,isol);
            z[isol] = m_pMat->GetSolute(isol)->ChargeNumber();
            M[isol] = m_pMat->GetSolute(isol)->MolarMass();
            d0[isol] = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
            kzc += z[isol]*spt.m_c[isol]*spt.m_k[isol];
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                d0p[isol][jsol] = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, jsol);
            }
        }
        
        //Get dk/dt (partial differential wrt time)
        vector<double> dkdt(nsol,0);
        for (int isol=0; isol<nsol; ++isol)
        {
            for (int jsol=0; jsol<nsol; ++jsol)
            {
                dkdt[isol] += spt.m_dkdc[isol][jsol]*spt.m_cdot[jsol];
            }
        }
        
        // evaluate the chat
        vector<double> vbardzdc(nsol, 0.0);
        vector<vector<double>> vdzdc(nsol, vector<double>(nsol, 0.0));
        
        // chemical reactions
        for (i=0; i<nreact; ++i) {
            double vbar = m_pMat->GetReaction(i)->m_Vbar;
            for (int isol = 0; isol < nsol; ++isol)
            {
                vbardzdc[isol] += vbar*m_pMat->GetReaction(i)->Tangent_ReactionSupply_Concentration(mp,isol);
                for (int jsol = 0; jsol < nsol; ++jsol)
                {
                    double v = m_pMat->GetReaction(i)->m_v[isol];
                    vdzdc[isol][jsol] = v*m_pMat->GetReaction(i)->Tangent_ReactionSupply_Concentration(mp,jsol);
                }
            }
        }
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i4=0; i<neln; ++i, i4 += ndpn)
        {
            for (j=0, j4 = 0; j<neln; ++j, j4 += ndpn)
            {
                mat3d Kvv = vdotTdotv(gradN[i], cv, gradN[j]);
                vec3d kJv = (pt.m_gradef*(H[i]/Jf) + gradN[i])*H[j];
                vec3d kvJ = (svJ*gradN[i])*H[j] + (gradN[j]*dep+pt.m_gradef*H[j]*d2ep)*H[i];
                double kJJ = (H[j]*(ksi/dt - dJoJ) + gradN[j]*pt.m_vft)*H[i]/Jf;
                
                ke[i4  ][j4  ] += Kvv(0,0)*detJ;
                ke[i4  ][j4+1] += Kvv(0,1)*detJ;
                ke[i4  ][j4+2] += Kvv(0,2)*detJ;
                ke[i4  ][j4+3] += kvJ.x*detJ;
                
                ke[i4+1][j4  ] += Kvv(1,0)*detJ;
                ke[i4+1][j4+1] += Kvv(1,1)*detJ;
                ke[i4+1][j4+2] += Kvv(1,2)*detJ;
                ke[i4+1][j4+3] += kvJ.y*detJ;
                
                ke[i4+2][j4  ] += Kvv(2,0)*detJ;
                ke[i4+2][j4+1] += Kvv(2,1)*detJ;
                ke[i4+2][j4+2] += Kvv(2,2)*detJ;
                ke[i4+2][j4+3] += kvJ.z*detJ;
                
                ke[i4+3][j4  ] += kJv.x*detJ;
                ke[i4+3][j4+1] += kJv.y*detJ;
                ke[i4+3][j4+2] += kJv.z*detJ;
                ke[i4+3][j4+3] += kJJ*detJ;
                
                for (int isol=0; isol<nsol; ++isol) {
                    vec3d kcv = gradN[i]*H[j]*(spt.m_k[isol]*spt.m_c[isol] + kzc);
                    vec3d kvc = vec3d(0);
                    double kJc = H[i]*H[j]*vbardzdc[isol]*spt.m_k[isol];
                    int irow = i4+4+isol;
                    
                    for(int jsol=0; jsol<nsol; ++jsol)
                    {
                        kJc += H[i]*H[j]*vbardzdc[isol]*spt.m_dkdc[jsol][isol]*spt.m_c[jsol];
                        double kcc = 0;
                        if (isol == jsol)
                        {
                            kvc += (spt.m_gradc[isol]*spt.m_dkdc[isol][isol]*H[j]+gradN[j]*spt.m_k[isol])*H[i]*R*T;
                            kcc += -(ksi/dt*spt.m_k[isol]+spt.m_dkdc[isol][isol]*spt.m_cdot[isol] + dkdt[isol] + spt.m_c[isol]*spt.m_dkdc[isol][isol]*ksi/dt)*H[i]*H[j] + vdzdc[isol][isol]*spt.m_k[isol]*H[i]*H[j] + gradN[i]*((-spt.m_gradc[isol]*(spt.m_dkdc[isol][isol]*d0[isol]+spt.m_k[isol]*d0p[isol][isol])+pt.m_vft*(spt.m_dkdc[isol][isol]*spt.m_c[isol]+spt.m_k[isol]))*H[j]-gradN[j]*spt.m_k[isol]*d0[isol]);
                            for (int ksol=0; ksol<nsol; ++ksol)
                            {
                                if(isol==ksol)
                                {
                                    kcc += H[i]*H[j]*vdzdc[isol][isol]*spt.m_dkdc[ksol][isol]*spt.m_c[ksol] + gradN[i]*((-spt.m_gradc[isol]*(spt.m_dkdc[isol][isol]*d0[isol]+spt.m_k[isol]*d0p[isol][isol])+pt.m_vft*(spt.m_dkdc[isol][isol]*spt.m_c[isol]+spt.m_k[isol]))*H[j]-gradN[j]*spt.m_k[isol]*d0[isol])*z[isol]*penalty;
                                }
                                else
                                {
                                    kcc += H[i]*H[j]*vdzdc[isol][isol]*spt.m_dkdc[ksol][isol]*spt.m_c[ksol] + gradN[i]*(-spt.m_gradc[ksol]*(spt.m_dkdc[ksol][isol]*d0[ksol]+spt.m_k[ksol]*d0p[ksol][isol])+pt.m_vft*spt.m_dkdc[ksol][isol]*spt.m_c[ksol])*H[j]*z[ksol]*penalty;
                                }
                            }
                        }
                        else
                        {
                            kvc += spt.m_gradc[jsol]*spt.m_dkdc[jsol][isol]*H[i]*H[j]*R*T;
                            kcc += -H[i]*H[j]*(spt.m_dkdc[isol][jsol]*spt.m_cdot[isol] + spt.m_c[isol]*spt.m_dkdc[isol][jsol]*ksi/dt) + H[i]*H[j]*vdzdc[isol][jsol]*spt.m_k[jsol] + gradN[i]*(-spt.m_gradc[isol]*(spt.m_dkdc[isol][jsol]*d0[isol]+spt.m_k[isol]*d0p[isol][jsol])+pt.m_vft*spt.m_dkdc[isol][jsol]*spt.m_c[isol])*H[j];
                            for (int ksol=0; ksol<nsol; ++ksol)
                            {
                                if(jsol==ksol)
                                {
                                    kcc += H[i]*H[j]*vdzdc[isol][jsol]*spt.m_dkdc[ksol][jsol]*spt.m_c[ksol] + gradN[i]*((-spt.m_gradc[jsol]*(spt.m_dkdc[jsol][jsol]*d0[jsol]+spt.m_k[jsol]*d0p[jsol][jsol])+pt.m_vft*(spt.m_dkdc[jsol][jsol]*spt.m_c[jsol]+spt.m_k[jsol]))*H[j]-gradN[j]*spt.m_k[jsol]*d0[jsol])*z[jsol]*penalty;
                                }
                                else
                                {
                                    kcc += H[i]*H[j]*vdzdc[isol][jsol]*spt.m_dkdc[ksol][jsol]*spt.m_c[ksol] + gradN[i]*(-spt.m_gradc[ksol]*(spt.m_dkdc[ksol][jsol]*d0[ksol]+spt.m_k[ksol]*d0p[ksol][jsol])+pt.m_vft*spt.m_dkdc[ksol][jsol]*spt.m_c[ksol])*H[j]*z[ksol]*penalty;
                                }
                            }
                        }
                        ke[irow][j4+4+jsol] += kcc*detJ;
                    }
                    ke[i4  ][j4+4+isol] += kvc.x*detJ*dms;
                    ke[i4+1][j4+4+isol] += kvc.y*detJ*dms;
                    ke[i4+2][j4+4+isol] += kvc.z*detJ*dms;
                    ke[i4+3][j4+4+isol] += kJc*detJ;
                    ke[irow][j4  ] += kcv.x*detJ;
                    ke[irow][j4+1] += kcv.y*detJ;
                    ke[irow][j4+2] += kcv.z*detJ;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::StiffnessMatrix(FELinearSystem& LS)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        // element stiffness matrix
        FEElementMatrix ke(el);
        
        // create the element's stiffness matrix
        int nsol = m_pMat->Solutes();
        int ndpn = 4 + nsol;
        int ndof = ndpn*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate material stiffness
        if (m_newFormulation) ElementStiffnessNew(el, ke);
        else ElementStiffness(el, ke);
        
        // get the element's LM vector
        vector<int> lm;
        UnpackLM(el, lm);
        ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
        LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::MassMatrix(FELinearSystem& LS)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        // element stiffness matrix
        FEElementMatrix ke(el);
        
        // create the element's stiffness matrix
        const int nsol = m_pMat->Solutes();
        const int ndpn = 4 + nsol;
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

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        // element stiffness matrix
        FEElementMatrix ke(el);
        
        // create the element's stiffness matrix
        const int nsol = m_pMat->Solutes();
        const int ndpn = 4 + nsol;
        int ndof = ndpn*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        if (m_newFormulation) ElementBodyForceStiffnessNew(bf, el, ke);
        else ElementBodyForceStiffness(bf, el, ke);
        
        // get the element's LM vector
        vector<int> lm;
        UnpackLM(el, lm);
        ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
        LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
//! calculates element inertial stiffness matrix
void FEFluidSolutesDomain3D::ElementMassMatrix(FESolidElement& el, matrix& ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i4, j, j4, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    const int ndpn = 4 + nsol;

    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double *H;
    double *Gr, *Gs, *Gt;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf)*m_btrans;
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjac0(el, Ji, n)*gw[n]*tp.alphaf;
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // setup the material point
        // NOTE: deformation gradient and determinant have already been evaluated in the stress routine
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        double dens = m_pMat->Fluid()->Density(mp);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i4=0; i<neln; ++i, i4 += ndpn)
        {
            for (j=0, j4 = 0; j<neln; ++j, j4 += ndpn)
            {
                mat3d Mv = ((mat3dd(ksi/dt) + pt.m_Lf)*H[j] + mat3dd(gradN[j]*pt.m_vft))*(H[i]*dens*detJ);
                vec3d mJ = pt.m_aft*(-H[i]*H[j]*dens/(pt.m_ef+1)*detJ);
                
                ke[i4  ][j4  ] += Mv(0,0);
                ke[i4  ][j4+1] += Mv(0,1);
                ke[i4  ][j4+2] += Mv(0,2);
                ke[i4  ][j4+3] += mJ.x;
                
                ke[i4+1][j4  ] += Mv(1,0);
                ke[i4+1][j4+1] += Mv(1,1);
                ke[i4+1][j4+2] += Mv(1,2);
                ke[i4+1][j4+3] += mJ.y;
                
                ke[i4+2][j4  ] += Mv(2,0);
                ke[i4+2][j4+1] += Mv(2,1);
                ke[i4+2][j4+2] += Mv(2,2);
                ke[i4+2][j4+3] += mJ.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::Update(const FETimeInfo& tp)
{
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
                // reset the logfile mode
                berr = true;
                if (NegativeJacobian::DoOutput()) feLogError(e.what());
            }
        }
    }
    
    if (berr) throw NegativeJacobianDetected();
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
void FEFluidSolutesDomain3D::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double alphaf = tp.alphaf;
    double alpham = tp.alpham;
    
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
    vec3d vt[NELN], vp[NELN];
    vec3d at[NELN], ap[NELN];
    double et[NELN], ep[NELN];
    double aet[NELN], aep[NELN];
    vector< vector<double> > ct(nsol, vector<double>(NELN));
    vector< vector<double> > cp(nsol, vector<double>(NELN));
    vector< vector<double> > act(nsol, vector<double>(NELN));
    vector< vector<double> > acp(nsol, vector<double>(NELN));
    vector<int> sid(nsol);
    for (int j=0; j<nsol; ++j) sid[j] = m_pMat->GetSolute(j)->GetSoluteDOF();
    for (int j=0; j<neln; ++j) {
        FENode& node = m_pMesh->Node(el.m_node[j]);
        vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
        vp[j] = node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
        at[j] = node.get_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
        ap[j] = node.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
        et[j] = node.get(m_dofEF);
        ep[j] = node.get_prev(m_dofEF);
        aet[j] = node.get(m_dofAEF);
        aep[j] = node.get_prev(m_dofAEF);
        for (int k=0; k<nsol; ++k) {
            ct[k][j] = node.get(m_dofC + sid[k]);
            cp[k][j] = node.get_prev(m_dofC + sid[k]);
            act[k][j] = node.get(m_dofAC + sid[k]);
            acp[k][j] = node.get_prev(m_dofAC + sid[k]);
        }
    }
    
    // loop over the integration points and update
    // velocity, velocity gradient, acceleration
    // stress and pressure at the integration point
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());

        // material point data
        pt.m_vft = el.Evaluate(vt, n)*alphaf + el.Evaluate(vp, n)*(1-alphaf);
        pt.m_Lf = gradient(el, vt, n)*alphaf + gradient(el, vp, n)*(1-alphaf);
        pt.m_aft = pt.m_Lf*pt.m_vft;
        if (m_btrans) pt.m_aft += el.Evaluate(at, n)*alpham + el.Evaluate(ap, n)*(1-alpham);
        pt.m_ef = el.Evaluate(et, n)*alphaf + el.Evaluate(ep, n)*(1-alphaf);
        pt.m_gradef = gradient(el, et, n)*alphaf + gradient(el, ep, n)*(1-alphaf);
        pt.m_efdot = pt.m_gradef*pt.m_vft;
        if (m_btrans) pt.m_efdot += el.Evaluate(aet, n)*alpham + el.Evaluate(aep, n)*(1-alpham);
        for (int isol=0; isol < nsol; ++isol) {
            spt.m_c[isol] = el.Evaluate(ct[isol], n)*alphaf + el.Evaluate(cp[isol], n)*(1-alphaf);
            spt.m_gradc[isol] = gradient(el, ct[isol], n)*alphaf + gradient(el, cp[isol], n)*(1-alphaf);
            spt.m_cdot[isol] = 0.0;
            if (m_btrans) spt.m_cdot[isol] += el.Evaluate(act[isol], n)*alpham + el.Evaluate(acp[isol], n)*(1-alpham);
        }
        
        m_pMat->PartitionCoefficientFunctions(mp, spt.m_k, spt.m_dkdJ, spt.m_dkdc);
        
        // calculate the stress at this material point
        pt.m_sf = m_pMat->Fluid()->Stress(mp);
        
        spt.m_pe = m_pMat->Fluid()->Pressure(mp);
        
        // calculate the solute flux and actual concentration
        for (int isol=0; isol < nsol; ++isol)
        {
            spt.m_j[isol] = m_newFormulation ? m_pMat->SoluteDiffusiveFlux(mp, isol) : m_pMat->SoluteFlux(mp, isol);
            spt.m_ca[isol] = m_pMat->ConcentrationActual(mp, isol);
        }
        
        // calculate the fluid pressure
        pt.m_pf = m_pMat->PressureActual(mp);
        
        spt.m_psi = m_pMat->ElectricPotential(mp);
        spt.m_Ie = m_pMat->CurrentDensity(mp);
        
        // update chemical reaction element data
        for (int j=0; j<m_pMat->Reactions(); ++j)
            m_pMat->GetReaction(j)->UpdateElementData(mp);
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::InertialForces(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        const int nsol = m_pMat->Solutes();
        const int ndpn = 4+nsol;
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

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::ElementInertialForce(FESolidElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian determinant
    double detJ;
    
    const double* H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;
    
    double*    gw = el.GaussWeights();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        double dens = m_pMat->Fluid()->Density(mp);
        
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n];
        
        H = el.H(n);
        
        for (i=0; i<neln; ++i)
        {
            vec3d f = pt.m_aft*(dens*H[i]);
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[ndpn*i  ] -= f.x*detJ;
            fe[ndpn*i+1] -= f.y*detJ;
            fe[ndpn*i+2] -= f.z*detJ;
        }
    }
}
