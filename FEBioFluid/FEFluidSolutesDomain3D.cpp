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
#include "FEFluidSolutesDomain3D.h"
#include "FEFluidSolver.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioFluidSolutes.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEFluidSolutesDomain3D::FEFluidSolutesDomain3D(FEModel* pfem) : FESolidDomain(pfem), FEFluidDomain(pfem), m_dofW(pfem), m_dofAW(pfem), m_dof(pfem)
{
    m_pMat = 0;
    m_btrans = true;
    
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
            ps.m_cdot.assign(nsol,0);
            ps.m_gradc.assign(nsol,vec3d(0,0,0));
            ps.m_j.assign(nsol,vec3d(0,0,0));
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
                ps.m_gradc[isol] = gradient(el, c0[isol], n);
            }
            
            for (int isol = 0; isol<nsol; ++isol)
                ps.m_j[isol] = m_pMat->SoluteFlux(mp, isol);
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
            
            if (pt.m_Jf <= 0) {
                feLogError("Negative jacobian was detected.");
                throw DoRunningRestart();
            }
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::InternalForces(FEGlobalVector& R, const FETimeInfo& tp)
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
        int nsol = m_pMat->Solutes();
        int ndpn = 4+nsol;
        int ndof = ndpn*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForce(el, fe, tp);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEFluidSolutesDomain3D::ElementInternalForce(FESolidElement& el, vector<double>& fe, const FETimeInfo& tp)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJ;
    
    mat3ds sv;
    vec3d gradp;
    
    const double *H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;

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
        gradp = pt.m_gradJf*m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        
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
        double dJoJ = pt.m_Jfdot/pt.m_Jf;
        
        for (i=0; i<neln; ++i)
        {
            vec3d fs = sv*gradN[i] + gradp*H[i];
            double fJ = dJoJ*H[i] + gradN[i]*pt.m_vft;
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[ndpn*i  ] -= fs.x*detJ;
            fe[ndpn*i+1] -= fs.y*detJ;
            fe[ndpn*i+2] -= fs.z*detJ;
            fe[ndpn*i+3] -= fJ*detJ;
            for (int isol=0; isol<nsol; ++isol)
                fe[ndpn*i+4+isol] -= (spt.m_j[isol]*gradN[i] - H[i]*(spt.m_cdot[isol]+spt.m_c[isol]*pt.m_Jfdot/pt.m_Jf))*detJ*dt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::BodyForce(FEGlobalVector& R, const FETimeInfo& tp, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int nsol = m_pMat->Solutes();
        int ndpn = 4+nsol;
        int ndof = ndpn*el.Nodes();
        fe.assign(ndof, 0);
        
        // apply body forces
        ElementBodyForce(BF, el, fe, tp);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEFluidSolutesDomain3D::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe, const FETimeInfo& tp)
{
    // jacobian
    double detJ;
    double *H;
    double* gw = el.GaussWeights();
    vec3d f;
    
    // number of nodes
    int neln = el.Nodes();
    int nsol = m_pMat->Solutes();
    int ndpn = 4+nsol;

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
        double dens = m_pMat->Fluid()->Density(mp);
        
        pt.m_r0 = el.Evaluate(r0, n);
        
        detJ = detJ0(el, n)*gw[n];
        
        // get the force
        f = BF.force(mp);
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i)
        {
            fe[ndpn*i  ] -= H[i]*dens*f.x*detJ;
            fe[ndpn*i+1] -= H[i]*dens*f.y*detJ;
            fe[ndpn*i+2] -= H[i]*dens*f.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEFluidSolutesDomain3D::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke, const FETimeInfo& tp)
{
    int neln = el.Nodes();
    int ndpn = ke.columns()/neln;
    
    // jacobian
    double detJ;
    double *H;
    double* gw = el.GaussWeights();
    vec3d f, k;
    
    // gradient of shape functions
    vec3d gradN;
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
        
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n];
        
        H = el.H(n);
        
        double dens = m_pMat->Fluid()->Density(mp);
        
        // get the force
        f = BF.force(mp);
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i) {
            for (int j=0; j<neln; ++j)
            {
                k = f*(-H[i]*H[j]*dens/pt.m_Jf*detJ);
                ke[ndpn*i  ][ndpn*j+3] += k.x;
                ke[ndpn*i+1][ndpn*j+3] += k.y;
                ke[ndpn*i+2][ndpn*j+3] += k.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEFluidSolutesDomain3D::ElementStiffness(FESolidElement &el, matrix &ke, const FETimeInfo& tp)
{
    int i, i4, j, j4, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    const int nsol = m_pMat->Solutes();
    const int ndpn = 4 + nsol;
    
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
        detJ = invjac0(el, Ji, n)*gw[n];
        
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

        // get the tangents
        mat3ds svJ = m_pMat->Fluid()->GetViscous()->Tangent_Strain(mp);
        tens4ds cv = m_pMat->Fluid()->Tangent_RateOfDeformation(mp);
        double dp = m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        double d2p = m_pMat->Fluid()->Tangent_Pressure_Strain_Strain(mp);
        // Jdot/J
        double dJoJ = pt.m_Jfdot/pt.m_Jf;
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i4=0; i<neln; ++i, i4 += ndpn)
        {
            for (j=0, j4 = 0; j<neln; ++j, j4 += ndpn)
            {
                mat3d Kvv = vdotTdotv(gradN[i], cv, gradN[j]);
                vec3d kJv = (pt.m_gradJf*(H[i]/pt.m_Jf) + gradN[i])*H[j];
                vec3d kvJ = (svJ*gradN[i])*H[j] + (gradN[j]*dp+pt.m_gradJf*(H[j]*d2p))*H[i];
                //                double kJJ = (H[j]*((ksi*m_btrans)/dt - dJoJ) + gradN[j]*pt.m_vft)*H[i]/pt.m_Jf;
                double kJJ = (H[j]*(ksi/dt - dJoJ) + gradN[j]*pt.m_vft)*H[i]/pt.m_Jf;
                
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
                    double d0 = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
                    double d0p = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, isol);
                    vec3d kcv = (spt.m_gradc[isol] + pt.m_gradJf*(spt.m_c[isol]/pt.m_Jf))*(-H[i]*H[j]);
                    double kcJ = -(H[i]*spt.m_c[isol]/pt.m_Jf*((ksi/dt-dJoJ)*H[j]+gradN[j]*pt.m_vft));
                    double kcc = -(H[i]*((ksi/dt+dJoJ)*H[j]+gradN[j]*pt.m_vft)+(gradN[j]*d0+spt.m_gradc[isol]*H[j]*d0p)*gradN[i]);
                    int irow = i4+4+isol;
                    ke[irow][j4  ] += kcv.x*detJ*dt;
                    ke[irow][j4+1] += kcv.y*detJ*dt;
                    ke[irow][j4+2] += kcv.z*detJ*dt;
                    ke[irow][j4+3] += kcJ*detJ*dt;
                    ke[irow][j4+4+isol] += kcc*detJ*dt;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
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
        ElementStiffness(el, ke, tp);
        
        // get the element's LM vector
        vector<int> lm;
        UnpackLM(el, lm);
        ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
        LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::MassMatrix(FELinearSystem& LS, const FETimeInfo& tp)
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
        ElementMassMatrix(el, ke, tp);
        
        // get the element's LM vector
        vector<int> lm;
        UnpackLM(el, lm);
        ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
        LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::BodyForceStiffness(FELinearSystem& LS, const FETimeInfo& tp, FEBodyForce& bf)
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
        ElementBodyForceStiffness(bf, el, ke, tp);
        
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
void FEFluidSolutesDomain3D::ElementMassMatrix(FESolidElement& el, matrix& ke, const FETimeInfo& tp)
{
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
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        
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
                mat3d Mv = ((mat3dd(ksi*m_btrans/dt) + pt.m_Lf)*H[j] + mat3dd(gradN[j]*pt.m_vft))*(H[i]*dens*detJ);
                vec3d mJ = pt.m_aft*(-H[i]*H[j]*dens/pt.m_Jf*detJ);
                
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
    
    // if we encountered an error, we request a running restart
    if (berr)
    {
        if (NegativeJacobian::DoOutput() == false) feLogError("Negative jacobian was detected.");
        throw DoRunningRestart();
    }
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
        pt.m_Jf = 1 + el.Evaluate(et, n)*alphaf + el.Evaluate(ep, n)*(1-alphaf);
        pt.m_gradJf = gradient(el, et, n)*alphaf + gradient(el, ep, n)*(1-alphaf);
        pt.m_Jfdot = pt.m_gradJf*pt.m_vft;
        if (m_btrans) pt.m_Jfdot += el.Evaluate(aet, n)*alpham + el.Evaluate(aep, n)*(1-alpham);
        for (int isol=0; isol < nsol; ++isol) {
            spt.m_c[isol] = el.Evaluate(ct[isol], n)*alphaf + el.Evaluate(cp[isol], n)*(1-alphaf);
            spt.m_gradc[isol] = gradient(el, ct[isol], n)*alphaf + gradient(el, cp[isol], n)*(1-alphaf);
            spt.m_cdot[isol] = spt.m_gradc[isol]*pt.m_vft;
            if (m_btrans) spt.m_cdot[isol] += el.Evaluate(act[isol], n)*alpham + el.Evaluate(acp[isol], n)*(1-alpham);
        }
        
        // calculate the stress at this material point
        pt.m_sf = m_pMat->Fluid()->Stress(mp);
        
        // calculate the fluid pressure
        pt.m_pf = m_pMat->Fluid()->Pressure(mp);
        
        // calculate the solute flux
        for (int isol=0; isol < nsol; ++isol)
            spt.m_j[isol] = m_pMat->SoluteFlux(mp, isol);
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::InertialForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    int NE = (int)m_Elem.size();
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
        ElementInertialForce(el, fe, tp);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain3D::ElementInertialForce(FESolidElement& el, vector<double>& fe, const FETimeInfo& tp)
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
