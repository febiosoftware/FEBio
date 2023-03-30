/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include "FEPolarFluidDomain3D.h"
#include "FEPolarFluidMaterialPoint.h"
#include "FEFluidSolver.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioPolarFluid.h"
#include <FECore/FELinearSystem.h>
#include "FEBodyMoment.h"

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEPolarFluidDomain3D::FEPolarFluidDomain3D(FEModel* pfem) : FESolidDomain(pfem), FEPolarFluidDomain(pfem), m_dofW(pfem), m_dofAW(pfem), m_dofG(pfem), m_dofAG(pfem), m_dof(pfem)
{
    m_pMat = 0;
    m_btrans = true;
    
    if (pfem)
    {
        m_dofW.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::RELATIVE_FLUID_VELOCITY));
        m_dofAW.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::RELATIVE_FLUID_ACCELERATION));
        
        m_dofG.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_ANGULAR_VELOCITY));
        m_dofAG.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_ANGULAR_ACCELERATION));
        
        m_dofEF = pfem->GetDOFIndex(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_DILATATION), 0);
        m_dofAEF = pfem->GetDOFIndex(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_DILATATION_TDERIV), 0);
        
        FEDofList dofs(pfem);
        dofs.AddDofs(m_dofW);
        dofs.AddDofs(m_dofG);
        dofs.AddDof(m_dofEF);
        m_dof = dofs;
    }
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEPolarFluidDomain3D& FEPolarFluidDomain3D::operator = (FEPolarFluidDomain3D& d)
{
    m_Elem = d.m_Elem;
    m_pMesh = d.m_pMesh;
    return (*this);
}

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::Serialize(DumpStream& ar)
{
    FESolidDomain::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_pMat;
    ar & m_dofW & m_dofAW;
    ar & m_dofG & m_dofAG;
    ar & m_dof;
    ar & m_dofEF & m_dofAEF;
}

//-----------------------------------------------------------------------------
// get total dof list
const FEDofList& FEPolarFluidDomain3D::GetDOFList() const
{
    return m_dof;
}

//-----------------------------------------------------------------------------
//! Assign material
void FEPolarFluidDomain3D::SetMaterial(FEMaterial* pmat)
{
    FEDomain::SetMaterial(pmat);
    if (pmat)
    {
        m_pMat = dynamic_cast<FEPolarFluid*>(pmat);
        assert(m_pMat);
    }
    else m_pMat = 0;
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEPolarFluidDomain3D::PreSolveUpdate(const FETimeInfo& timeInfo)
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
            
            if (pt.m_ef <= -1) {
                throw NegativeJacobianDetected();
            }
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::InternalForces(FEGlobalVector& R)
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
        int ndof = 7*el.Nodes();
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

void FEPolarFluidDomain3D::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJ;
    
    const double *H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double*    gw = el.GaussWeights();
    
    FEViscousFluid* Viscous = m_pMat->GetViscous();
    FEViscousPolarFluid* ViscPol = m_pMat->GetViscousPolar();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        // get the viscous stress tensor for this integration point
        vec3d th = ViscPol->SkewStressDualVector(mp);
        mat3d sva = mat3da(th);
        mat3d sv = sva + Viscous->Stress(mp);
        mat3d M = ViscPol->CoupleStress(mp);
        // get the gradient of the elastic pressure
        vec3d gradp = pt.m_gradef*m_pMat->Tangent_Pressure_Strain(mp);
        
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
        
        for (i=0; i<neln; ++i)
        {
            vec3d fv = sv*gradN[i] + gradp*H[i];
            vec3d fg = M*gradN[i] - th*(2*H[i]);
            double fJ = dJoJ*H[i] + gradN[i]*pt.m_vft;
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[7*i  ] -= fv.x*detJ;
            fe[7*i+1] -= fv.y*detJ;
            fe[7*i+2] -= fv.z*detJ;
            fe[7*i+3] -= fg.x*detJ;
            fe[7*i+4] -= fg.y*detJ;
            fe[7*i+5] -= fg.z*detJ;
            fe[7*i+6] -= fJ*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 7*el.Nodes();
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

void FEPolarFluidDomain3D::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
{
    // jacobian
    double detJ;
    double *H;
    double* gw = el.GaussWeights();
    vec3d f;
    
    // number of nodes
    int neln = el.Nodes();
    
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
        double dens = m_pMat->Density(mp);
        
        pt.m_r0 = el.Evaluate(r0, n);
        
        detJ = detJ0(el, n)*gw[n];
        
        // get the force
        f = BF.force(mp);
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i)
        {
            fe[7*i  ] -= H[i]*dens*f.x*detJ;
            fe[7*i+1] -= H[i]*dens*f.y*detJ;
            fe[7*i+2] -= H[i]*dens*f.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEPolarFluidDomain3D::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int neln = el.Nodes();
    
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
        detJ = detJ0(el, n)*gw[n]*tp.alphaf;
        
        H = el.H(n);
        
        double dens = m_pMat->Density(mp);
        
        // get the force
        f = BF.force(mp);
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i) {
            for (int j=0; j<neln; ++j)
            {
                k = f*(-H[i]*H[j]*dens/(pt.m_ef+1)*detJ);
                ke[7*i  ][7*j+6] += k.x;
                ke[7*i+1][7*j+6] += k.y;
                ke[7*i+2][7*j+6] += k.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::BodyMoment(FEGlobalVector& R, FEBodyMoment& bm)
{
    int NE = (int)m_Elem.size();
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 7*el.Nodes();
        fe.assign(ndof, 0);
        
        // apply body forces
        ElementBodyMoment(bm, el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEPolarFluidDomain3D::ElementBodyMoment(FEBodyMoment& bm, FESolidElement& el, vector<double>& fe)
{
    // jacobian
    double detJ;
    double *H;
    double* gw = el.GaussWeights();
    vec3d m;
    
    // number of nodes
    int neln = el.Nodes();
    
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
        double dens = m_pMat->Density(mp);
        
        pt.m_r0 = el.Evaluate(r0, n);
        
        detJ = detJ0(el, n)*gw[n];
        
        // get the moment
        m = bm.moment(mp);
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i)
        {
            fe[7*i+3] -= H[i]*dens*m.x*detJ;
            fe[7*i+4] -= H[i]*dens*m.y*detJ;
            fe[7*i+5] -= H[i]*dens*m.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body moments
void FEPolarFluidDomain3D::ElementBodyMomentStiffness(FEBodyMoment& bm, FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int neln = el.Nodes();
    
    // jacobian
    double detJ;
    double *H;
    double* gw = el.GaussWeights();
    vec3d m, k;
    
    // gradient of shape functions
    vec3d gradN;
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
        
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n]*tp.alphaf;
        
        H = el.H(n);
        
        double dens = m_pMat->Density(mp);
        
        // get the force
        m = bm.moment(mp);
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i) {
            for (int j=0; j<neln; ++j)
            {
                k = m*(-H[i]*H[j]*dens/(pt.m_ef+1)*detJ);
                ke[7*i+3][7*j+6] += k.x;
                ke[7*i+4][7*j+6] += k.y;
                ke[7*i+5][7*j+6] += k.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEPolarFluidDomain3D::ElementStiffness(FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i7, j, j7, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);
    
    double *H, *Gr, *Gs, *Gt;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    FEViscousFluid* Viscous = m_pMat->GetViscous();
    FEViscousPolarFluid* ViscPol = m_pMat->GetViscousPolar();
    
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
        double Jf = 1 + pt.m_ef;
        
        // get the tangents
        mat3d svJ = Viscous->Tangent_Strain(mp) + ViscPol->SkewTangent_Strain(mp);
        vec3d thJ = ViscPol->SkewTangent_Strain(mp).vec();
        tens4ds cvs = m_pMat->Tangent_RateOfDeformation(mp);
        mat3d cva = ViscPol->SkewTangent_RateOfRotation(mp);
        tens4d m = ViscPol->CoupleTangent_RateOfRotation(mp);
        mat3d mJ = ViscPol->CoupleTangent_Strain(mp);
        double dp = m_pMat->Tangent_Pressure_Strain(mp);
        double d2p = m_pMat->Tangent_Pressure_Strain_Strain(mp);
        // Jdot/J
        double dJoJ = pt.m_efdot/Jf;
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i7=0; i<neln; ++i, i7 += 7)
        {
            for (j=0, j7 = 0; j<neln; ++j, j7 += 7)
            {
                mat3da gNa(gradN[i]), gNb(gradN[j]);
                mat3d Kvv = vdotTdotv(gradN[i], cvs, gradN[j]) + gNa*cva*gNb/2;
                mat3d Kgv = cva*gNb*H[i];
                vec3d kJv = (pt.m_gradef*(H[i]/Jf) + gradN[i])*H[j];
                mat3d Kvg = -gNa*cva*H[j];
                //! TODO: There may be a mistake in vdotTdotv when using tens4d object, evaluates as b.T.a instead of a.T.b
//                mat3d Kgg = vdotTdotv(gradN[i], m, gradN[j])-cva*(2*H[i]*H[j]);
                mat3d Kgg = vdotTdotv(gradN[j], m, gradN[i])-cva*(2*H[i]*H[j]);
                vec3d kvJ = (svJ*gradN[i])*H[j] + (gradN[j]*dp+pt.m_gradef*(H[j]*d2p))*H[i];
                vec3d kgJ = mJ*gradN[i]*H[j] - thJ*(2*H[i]*H[j]);
                double kJJ = (H[j]*(ksi/dt - dJoJ) + gradN[j]*pt.m_vft)*H[i]/Jf;
                
                ke[i7  ][j7  ] += Kvv(0,0)*detJ;
                ke[i7  ][j7+1] += Kvv(0,1)*detJ;
                ke[i7  ][j7+2] += Kvv(0,2)*detJ;
                ke[i7  ][j7+3] += Kvg(0,0)*detJ;
                ke[i7  ][j7+4] += Kvg(0,1)*detJ;
                ke[i7  ][j7+5] += Kvg(0,2)*detJ;
                ke[i7  ][j7+6] += kvJ.x*detJ;
                
                ke[i7+1][j7  ] += Kvv(1,0)*detJ;
                ke[i7+1][j7+1] += Kvv(1,1)*detJ;
                ke[i7+1][j7+2] += Kvv(1,2)*detJ;
                ke[i7+1][j7+3] += Kvg(1,0)*detJ;
                ke[i7+1][j7+4] += Kvg(1,1)*detJ;
                ke[i7+1][j7+5] += Kvg(1,2)*detJ;
                ke[i7+1][j7+6] += kvJ.y*detJ;
                
                ke[i7+2][j7  ] += Kvv(2,0)*detJ;
                ke[i7+2][j7+1] += Kvv(2,1)*detJ;
                ke[i7+2][j7+2] += Kvv(2,2)*detJ;
                ke[i7+2][j7+3] += Kvg(2,0)*detJ;
                ke[i7+2][j7+4] += Kvg(2,1)*detJ;
                ke[i7+2][j7+5] += Kvg(2,2)*detJ;
                ke[i7+2][j7+6] += kvJ.z*detJ;
                
                ke[i7+3][j7  ] += Kgv(0,0)*detJ;
                ke[i7+3][j7+1] += Kgv(0,1)*detJ;
                ke[i7+3][j7+2] += Kgv(0,2)*detJ;
                ke[i7+3][j7+3] += Kgg(0,0)*detJ;
                ke[i7+3][j7+4] += Kgg(0,1)*detJ;
                ke[i7+3][j7+5] += Kgg(0,2)*detJ;
                ke[i7+3][j7+6] += kgJ.x*detJ;
                
                ke[i7+4][j7  ] += Kgv(1,0)*detJ;
                ke[i7+4][j7+1] += Kgv(1,1)*detJ;
                ke[i7+4][j7+2] += Kgv(1,2)*detJ;
                ke[i7+4][j7+3] += Kgg(1,0)*detJ;
                ke[i7+4][j7+4] += Kgg(1,1)*detJ;
                ke[i7+4][j7+5] += Kgg(1,2)*detJ;
                ke[i7+4][j7+6] += kgJ.y*detJ;
                
                ke[i7+5][j7  ] += Kgv(2,0)*detJ;
                ke[i7+5][j7+1] += Kgv(2,1)*detJ;
                ke[i7+5][j7+2] += Kgv(2,2)*detJ;
                ke[i7+5][j7+3] += Kgg(2,0)*detJ;
                ke[i7+5][j7+4] += Kgg(2,1)*detJ;
                ke[i7+5][j7+5] += Kgg(2,2)*detJ;
                ke[i7+5][j7+6] += kgJ.z*detJ;
                
                ke[i7+6][j7  ] += kJv.x*detJ;
                ke[i7+6][j7+1] += kJv.y*detJ;
                ke[i7+6][j7+2] += kJv.z*detJ;
                ke[i7+6][j7+6] += kJJ*detJ;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::StiffnessMatrix(FELinearSystem& LS)
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
        int ndof = 7*el.Nodes();
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

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::MassMatrix(FELinearSystem& LS)
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
        int ndof = 7*el.Nodes();
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
void FEPolarFluidDomain3D::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
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
        int ndof = 7*el.Nodes();
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

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::BodyMomentStiffness(FELinearSystem& LS, FEBodyMoment& bm)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        // element stiffness matrix
        FEElementMatrix ke(el);
        
        // create the element's stiffness matrix
        int ndof = 7*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementBodyMomentStiffness(bm, el, ke);
        
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
void FEPolarFluidDomain3D::ElementMassMatrix(FESolidElement& el, matrix& ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i7, j, j7, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
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
        FEPolarFluidMaterialPoint& pf = *(mp.ExtractData<FEPolarFluidMaterialPoint>());
        
        double dens = m_pMat->Density(mp);
        double kg = m_pMat->m_kg;
        double moi = dens*kg*kg;
        double J = 1+pt.m_ef;
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i7=0; i<neln; ++i, i7 += 7)
        {
            for (j=0, j7 = 0; j<neln; ++j, j7 += 7)
            {
                mat3d Mv = ((mat3dd(ksi/dt) + pt.m_Lf)*H[j] + mat3dd(gradN[j]*pt.m_vft))*(H[i]*dens*detJ);
                mat3d Jgg = mat3dd(ksi/dt*H[j] + gradN[j]*pt.m_vft)*(H[i]*moi*detJ);
                mat3d Jgv = pf.m_Psi*(H[i]*H[j]*moi*detJ);
                vec3d mJ = pt.m_aft*(-H[i]*H[j]*dens/J*detJ);
                vec3d jJ = pf.m_gfdot*(-H[i]*H[j]*moi/J*detJ);
                
                ke[i7  ][j7  ] += Mv(0,0);
                ke[i7  ][j7+1] += Mv(0,1);
                ke[i7  ][j7+2] += Mv(0,2);
                ke[i7  ][j7+6] += mJ.x;
                
                ke[i7+1][j7  ] += Mv(1,0);
                ke[i7+1][j7+1] += Mv(1,1);
                ke[i7+1][j7+2] += Mv(1,2);
                ke[i7+1][j7+6] += mJ.y;
                
                ke[i7+2][j7  ] += Mv(2,0);
                ke[i7+2][j7+1] += Mv(2,1);
                ke[i7+2][j7+2] += Mv(2,2);
                ke[i7+2][j7+6] += mJ.z;

                ke[i7+3][j7  ] += Jgv(0,0);
                ke[i7+3][j7+1] += Jgv(0,1);
                ke[i7+3][j7+2] += Jgv(0,2);
                ke[i7+3][j7+3] += Jgg(0,0);
                ke[i7+3][j7+4] += Jgg(0,1);
                ke[i7+3][j7+5] += Jgg(0,2);
                ke[i7+3][j7+6] += jJ.x;
                
                ke[i7+4][j7  ] += Jgv(1,0);
                ke[i7+4][j7+1] += Jgv(1,1);
                ke[i7+4][j7+2] += Jgv(1,2);
                ke[i7+4][j7+3] += Jgg(1,0);
                ke[i7+4][j7+4] += Jgg(1,1);
                ke[i7+4][j7+5] += Jgg(1,2);
                ke[i7+4][j7+6] += jJ.y;
                
                ke[i7+5][j7  ] += Jgv(2,0);
                ke[i7+5][j7+1] += Jgv(2,1);
                ke[i7+5][j7+2] += Jgv(2,2);
                ke[i7+5][j7+3] += Jgg(2,0);
                ke[i7+5][j7+4] += Jgg(2,1);
                ke[i7+5][j7+5] += Jgg(2,2);
                ke[i7+5][j7+6] += jJ.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::Update(const FETimeInfo& tp)
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
void FEPolarFluidDomain3D::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double alphaf = tp.alphaf;
    double alpham = tp.alpham;
    
    // get the solid element
    FESolidElement& el = m_Elem[iel];
    
    // get the number of integration points
    int nint = el.GaussPoints();
    
    // number of nodes
    int neln = el.Nodes();
    
    // nodal coordinates
    const int NELN = FEElement::MAX_NODES;
    vec3d vt[NELN], vp[NELN];
    vec3d at[NELN], ap[NELN];
    vec3d gt[NELN], gp[NELN];
    vec3d agt[NELN], agp[NELN];
    double et[NELN], ep[NELN];
    double aet[NELN], aep[NELN];
    for (int j=0; j<neln; ++j) {
        FENode& node = m_pMesh->Node(el.m_node[j]);
        vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
        vp[j] = node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
        at[j] = node.get_vec3d(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
        ap[j] = node.get_vec3d_prev(m_dofAW[0], m_dofAW[1], m_dofAW[2]);
        gt[j] = node.get_vec3d(m_dofG[0], m_dofG[1], m_dofG[2]);
        gp[j] = node.get_vec3d_prev(m_dofG[0], m_dofG[1], m_dofG[2]);
        agt[j] = node.get_vec3d(m_dofAG[0], m_dofAG[1], m_dofAG[2]);
        agp[j] = node.get_vec3d_prev(m_dofAG[0], m_dofAG[1], m_dofAG[2]);
        et[j] = node.get(m_dofEF);
        ep[j] = node.get_prev(m_dofEF);
        aet[j] = node.get(m_dofAEF);
        aep[j] = node.get_prev(m_dofAEF);
    }
    
    // loop over the integration points and update
    // velocity, velocity gradient, acceleration
    // stress and pressure at the integration point
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEPolarFluidMaterialPoint& pf = *(mp.ExtractData<FEPolarFluidMaterialPoint>());

        // material point data
        pt.m_vft = el.Evaluate(vt, n)*alphaf + el.Evaluate(vp, n)*(1-alphaf);
        pt.m_Lf = gradient(el, vt, n)*alphaf + gradient(el, vp, n)*(1-alphaf);
        pt.m_aft = pt.m_Lf*pt.m_vft;
        if (m_btrans) pt.m_aft += el.Evaluate(at, n)*alpham + el.Evaluate(ap, n)*(1-alpham);
        pf.m_gf = el.Evaluate(gt, n)*alphaf + el.Evaluate(gp, n)*(1-alphaf);
        pf.m_Psi = gradient(el, gt, n)*alphaf + gradient(el, gp, n)*(1-alphaf);
        pf.m_gfdot = pf.m_Psi*pt.m_vft;
        if (m_btrans) pf.m_gfdot += el.Evaluate(agt, n)*alpham + el.Evaluate(agp, n)*(1-alpham);
        pt.m_ef = el.Evaluate(et, n)*alphaf + el.Evaluate(ep, n)*(1-alphaf);
        pt.m_gradef = gradient(el, et, n)*alphaf + gradient(el, ep, n)*(1-alphaf);
        pt.m_efdot = pt.m_gradef*pt.m_vft;
        if (m_btrans) pt.m_efdot += el.Evaluate(aet, n)*alpham + el.Evaluate(aep, n)*(1-alpham);
        
        // calculate the stress at this material point
        pt.m_sf = m_pMat->Stress(mp);
        pf.m_sfa = m_pMat->GetViscousPolar()->SkewStress(mp);
        
        // calculate the fluid pressure
        pt.m_pf = m_pMat->Pressure(mp);
    }
}

//-----------------------------------------------------------------------------
void FEPolarFluidDomain3D::InertialForces(FEGlobalVector& R)
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
        int ndof = 7*el.Nodes();
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
void FEPolarFluidDomain3D::ElementInertialForce(FESolidElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian determinant
    double detJ;
    
    const double* H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*    gw = el.GaussWeights();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEPolarFluidMaterialPoint& pf = *(mp.ExtractData<FEPolarFluidMaterialPoint>());
        double dens = m_pMat->Density(mp);
        double kg = m_pMat->m_kg;
        
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n];
        
        H = el.H(n);
        
        for (i=0; i<neln; ++i)
        {
            vec3d fv = pt.m_aft*(dens*H[i]);
            vec3d fg = pf.m_gfdot*(dens*kg*kg*H[i]);
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[7*i  ] -= fv.x*detJ;
            fe[7*i+1] -= fv.y*detJ;
            fe[7*i+2] -= fv.z*detJ;
            fe[7*i+3] -= fg.x*detJ;
            fe[7*i+4] -= fg.y*detJ;
            fe[7*i+5] -= fg.z*detJ;
        }
    }
}
