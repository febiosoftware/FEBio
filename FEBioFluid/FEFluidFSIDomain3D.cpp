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
#include "FEFluidFSIDomain3D.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include "FEBioFSI.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEFluidFSIDomain3D::FEFluidFSIDomain3D(FEModel* pfem) : FESolidDomain(pfem), FEFluidFSIDomain(pfem), m_dofU(pfem), m_dofV(pfem), m_dofW(pfem), m_dofAW(pfem), m_dofSU(pfem), m_dofR(pfem), m_dof(pfem)
{
    m_pMat = 0;
    m_btrans = true;
    m_sseps = 0;

	if (pfem)
	{
		m_dofU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT));
		m_dofV.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::VELOCITY));
		m_dofW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY));
		m_dofAW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_ACCELERATION));
		m_dofSU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_DISPLACEMENT));
		m_dofR.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RIGID_ROTATION));
		m_dofEF = pfem->GetDOFIndex(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION), 0);
		m_dofAEF = pfem->GetDOFIndex(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION_TDERIV), 0);
	}
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEFluidFSIDomain3D& FEFluidFSIDomain3D::operator = (FEFluidFSIDomain3D& d)
{
    m_Elem = d.m_Elem;
    m_pMesh = d.m_pMesh;
    return (*this);
}

//-----------------------------------------------------------------------------
// get the total dof
const FEDofList& FEFluidFSIDomain3D::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
//! Assign material
void FEFluidFSIDomain3D::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
    if (pmat)
    {
        m_pMat = dynamic_cast<FEFluidFSI*>(pmat);
        assert(m_pMat);
    }
    else m_pMat = 0;
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::Activate()
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
            node.set_active(m_dofW[0]);
            node.set_active(m_dofW[1]);
            node.set_active(m_dofW[2]);
            node.set_active(m_dofEF);
        }
    }
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEFluidFSIDomain3D::PreSolveUpdate(const FETimeInfo& timeInfo)
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
                FEFSIMaterialPoint& ft = *mp.ExtractData<FEFSIMaterialPoint>();
                et.m_Wp = et.m_Wt;
                
                if ((pt.m_ef <= -1) || (et.m_J <= 0)) {
                    throw NegativeJacobianDetected();
                }
                
                mp.Update(timeInfo);
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data.
void FEFluidFSIDomain3D::UnpackLM(FEElement& el, vector<int>& lm)
{
    int N = el.Nodes();
    lm.resize(N*10);
    for (int i=0; i<N; ++i)
    {
        FENode& node = m_pMesh->Node(el.m_node[i]);
        vector<int>& id = node.m_ID;
        
        // first the displacement dofs
        lm[7*i  ] = id[m_dofU[0]];
        lm[7*i+1] = id[m_dofU[1]];
        lm[7*i+2] = id[m_dofU[2]];
        lm[7*i+3] = id[m_dofW[0]];
        lm[7*i+4] = id[m_dofW[1]];
        lm[7*i+5] = id[m_dofW[2]];
        lm[7*i+6] = id[m_dofEF];

		// rigid rotational dofs
		lm[7*N + 3*i  ] = id[m_dofR[0]];
		lm[7*N + 3*i+1] = id[m_dofR[1]];
		lm[7*N + 3*i+2] = id[m_dofR[2]];
    }
    
    // substitute interface dofs for solid-shell interfaces
    FESolidElement& sel = static_cast<FESolidElement&>(el);
    for (int i = 0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = m_pMesh->Node(el.m_node[i]);
            vector<int>& id = node.m_ID;
            
            // first the displacement dofs
            lm[7*i  ] = id[m_dofSU[0]];
            lm[7*i+1] = id[m_dofSU[1]];
            lm[7*i+2] = id[m_dofSU[2]];
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::InternalForces(FEGlobalVector& R)
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
        
        if (el.isActive()) {
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
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEFluidFSIDomain3D::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJ;
    
    mat3ds sv, ss;
    vec3d gradp;
    
    const double *H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double*	gw = el.GaussWeights();
    
    double dtrans = m_btrans ? 1 : m_sseps;
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEElasticMaterialPoint& et = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEFSIMaterialPoint& ft = *(mp.ExtractData<FEFSIMaterialPoint>());
        double Jf = 1 + pt.m_ef;
        
        // calculate the jacobian
        detJ = invjact(el, Ji, n, tp.alphaf)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        // get the viscous stress tensor for this integration point
        sv = m_pMat->Fluid()->GetViscous()->Stress(mp);
        // get the gradient of the elastic pressure
        gradp = pt.m_gradef*m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        // get the solid stress tensor
        ss = ft.m_ss;

        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
        {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        // Jfdot/Jf
        double dJfoJ = pt.m_efdot/Jf;
        // Jsdot/Js
        double dJsoJ = ft.m_Jdot/et.m_J;
        
        for (i=0; i<neln; ++i)
        {
            vec3d fs = (ss*gradN[i])*detJ;
            vec3d ff = (sv*gradN[i] + gradp*H[i])*detJ;
            double fJ = (H[i]*(((dJfoJ - dJsoJ)*dtrans + (pt.m_gradef*ft.m_w)/Jf)) + gradN[i]*ft.m_w)*detJ;
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[7*i  ] -= fs.x;
            fe[7*i+1] -= fs.y;
            fe[7*i+2] -= fs.z;
            fe[7*i+3] -= ff.x;
            fe[7*i+4] -= ff.y;
            fe[7*i+5] -= ff.z;
            fe[7*i+6] -= fJ;
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
    for (int i=0; i<NE; ++i)
    {
        // get the element
        FESolidElement& el = m_Elem[i];
        
        if (el.isActive()) {
            vector<double> fe;
            vector<int> lm;
            
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
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEFluidFSIDomain3D::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    // jacobian
    double detJ;
    double *H;
    double* gw = el.GaussWeights();
    vec3d f;
    
    // number of nodes
    int neln = el.Nodes();
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        double dens = m_pMat->Fluid()->Density(mp);
        
        detJ = detJt(el, n, tp.alphaf)*gw[n];
        
        // get the force
        f = BF.force(mp)*(dens*detJ);
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i)
        {
            fe[7*i+3] -= H[i]*f.x;
            fe[7*i+4] -= H[i]*f.y;
            fe[7*i+5] -= H[i]*f.z;
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
//! For now, we assume that the body force is constant
void FEFluidFSIDomain3D::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int neln = el.Nodes();
    int ndof = ke.columns()/neln;
    
    // jacobian
    double Ji[3][3], detJ;
    double *H, *Gr, *Gs, *Gt;
    double* gw = el.GaussWeights();
    vec3d f, kwJ;
    mat3d Kwu;
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
        
        // calculate the jacobian
        detJ = invjact(el, Ji, n, tp.alphaf)*gw[n]*tp.alphaf;
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        
        double dens = m_pMat->Fluid()->Density(mp);
        
        // get the force
        f = BF.force(mp)*(dens*detJ);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (int i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        for (int i=0; i<neln; ++i) {
            for (int j=0; j<neln; ++j)
            {
                kwJ = f*(-H[i]*H[j]/(pt.m_ef+1));
                Kwu = (f & gradN[j])*H[i];
                ke[ndof*i+3][ndof*j  ] += Kwu(0,0); ke[ndof*i+3][ndof*j+1] += Kwu(0,1); ke[ndof*i+3][ndof*j+2] += Kwu(0,2);
                ke[ndof*i+4][ndof*j  ] += Kwu(1,0); ke[ndof*i+4][ndof*j+1] += Kwu(1,1); ke[ndof*i+4][ndof*j+2] += Kwu(1,2);
                ke[ndof*i+5][ndof*j  ] += Kwu(2,0); ke[ndof*i+5][ndof*j+1] += Kwu(2,1); ke[ndof*i+5][ndof*j+2] += Kwu(2,2);
                ke[ndof*i+3][ndof*j+6] += kwJ.x;
                ke[ndof*i+4][ndof*j+6] += kwJ.y;
                ke[ndof*i+5][ndof*j+6] += kwJ.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEFluidFSIDomain3D::ElementStiffness(FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i7, j, j7, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double dt = tp.timeIncrement;
    double a = tp.gamma/(tp.beta*dt);
    double b = tp.alpham/(tp.alphaf*tp.beta*dt*dt);
    double c = tp.alpham/(tp.alphaf*tp.gamma*dt);

    double dtrans = m_btrans ? 1 : m_sseps;
    
    double *H, *Gr, *Gs, *Gt;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjact(el, Ji, n, tp.alphaf)*gw[n]*tp.alphaf;
        
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
        FEElasticMaterialPoint& et = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFSIMaterialPoint& fpt = *(mp.ExtractData<FEFSIMaterialPoint>());
        double Jf = 1 + pt.m_ef;
        
        // get the tangents
        mat3ds ss = m_pMat->Solid()->Stress(mp);
        tens4ds cs = m_pMat->Solid()->Tangent(mp);
        mat3ds sv = m_pMat->Fluid()->GetViscous()->Stress(mp);
        mat3ds svJ = m_pMat->Fluid()->GetViscous()->Tangent_Strain(mp);
        tens4ds cv = m_pMat->Fluid()->Tangent_RateOfDeformation(mp);
        double dp = m_pMat->Fluid()->Tangent_Pressure_Strain(mp);
        double d2p = m_pMat->Fluid()->Tangent_Pressure_Strain_Strain(mp);
        // Jfdot/Jf
        double dJfoJ = pt.m_efdot/Jf;
        vec3d gradp = pt.m_gradef*dp;
        // Jsdot/Js = div(vs)
        double dJsoJ = fpt.m_Jdot/et.m_J;

        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i7=0; i<neln; ++i, i7 += 7)
        {
            for (j=0, j7 = 0; j<neln; ++j, j7 += 7)
            {
                mat3d M = mat3dd(a*dtrans) - et.m_L;
                mat3d Kuu = (vdotTdotv(gradN[i], cs, gradN[j])
                             + mat3dd(gradN[i]*(ss*gradN[j])))*detJ;
                mat3d Kwu = (vdotTdotv(gradN[i], cv, gradN[j])*M
                             + sv*((gradN[i] & gradN[j]) - (gradN[j] & gradN[i]))
                             + ((gradp & gradN[j]) - (gradN[j] & gradp)*H[i]))*detJ;
                mat3d Kww = vdotTdotv(gradN[i], cv, gradN[j])*detJ;
                vec3d kwJ = ((svJ*gradN[i])*H[j]
                             + (gradN[j]*dp+pt.m_gradef*(H[j]*d2p))*H[i])*detJ;
                vec3d kJu = ((gradN[j]*(dJfoJ - dJsoJ - a)*dtrans
                              + ((gradN[j] & pt.m_gradef) - (pt.m_gradef & gradN[j]))*fpt.m_w/Jf + et.m_L.transpose()*gradN[j])*H[i]
                             + ((gradN[j] & gradN[i]) - (gradN[i] & gradN[j]))*fpt.m_w
                             )*detJ;
                vec3d kJw = ((pt.m_gradef*(H[i]/Jf) + gradN[i])*H[j])*detJ;
                double kJJ = (((c - dJfoJ)*dtrans - (pt.m_gradef*fpt.m_w)/Jf)*H[j]
                              + gradN[j]*fpt.m_w)*H[i]/Jf*detJ;
                
                ke[i7  ][j7  ] += Kuu(0,0); ke[i7  ][j7+1] += Kuu(0,1); ke[i7  ][j7+2] += Kuu(0,2);
                ke[i7+1][j7  ] += Kuu(1,0); ke[i7+1][j7+1] += Kuu(1,1); ke[i7+1][j7+2] += Kuu(1,2);
                ke[i7+2][j7  ] += Kuu(2,0); ke[i7+2][j7+1] += Kuu(2,1); ke[i7+2][j7+2] += Kuu(2,2);
                
                ke[i7+3][j7  ] += Kwu(0,0); ke[i7+3][j7+1] += Kwu(0,1); ke[i7+3][j7+2] += Kwu(0,2);
                ke[i7+4][j7  ] += Kwu(1,0); ke[i7+4][j7+1] += Kwu(1,1); ke[i7+4][j7+2] += Kwu(1,2);
                ke[i7+5][j7  ] += Kwu(2,0); ke[i7+5][j7+1] += Kwu(2,1); ke[i7+5][j7+2] += Kwu(2,2);
                
                ke[i7+3][j7+3] += Kww(0,0); ke[i7+3][j7+4] += Kww(0,1); ke[i7+3][j7+5] += Kww(0,2);
                ke[i7+4][j7+3] += Kww(1,0); ke[i7+4][j7+4] += Kww(1,1); ke[i7+4][j7+5] += Kww(1,2);
                ke[i7+5][j7+3] += Kww(2,0); ke[i7+5][j7+4] += Kww(2,1); ke[i7+5][j7+5] += Kww(2,2);
                
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
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::StiffnessMatrix(FELinearSystem& LS)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        if (el.isActive()) {
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
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::MassMatrix(FELinearSystem& LS)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        if (el.isActive()) {

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
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
    FEFluidFSI* pme = dynamic_cast<FEFluidFSI*>(GetMaterial()); assert(pme);
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        if (el.isActive()) {

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
}

//-----------------------------------------------------------------------------
//! calculates element inertial stiffness matrix
void FEFluidFSIDomain3D::ElementMassMatrix(FESolidElement& el, matrix& ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i7, j, j7, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    double dtrans = m_btrans ? 1 : m_sseps;

    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double *H;
    double *Gr, *Gs, *Gt;
    
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
        FEFSIMaterialPoint& fpt = *(mp.ExtractData<FEFSIMaterialPoint>());
        
        double dens = m_pMat->Fluid()->Density(mp);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i7=0; i<neln; ++i, i7 += 7)
        {
            for (j=0, j7 = 0; j<neln; ++j, j7 += 7)
            {
                mat3d Kwu = ((pt.m_aft & gradN[j]) - pt.m_Lf*(fpt.m_w*gradN[j])
                             + mat3dd((gradN[j]*fpt.m_w)*a*dtrans)
                             + mat3dd(b*H[j]*dtrans)
                             )*(H[i]*dens*detJ);
                mat3d Kww = (mat3dd(c*dtrans*H[j] + gradN[j]*fpt.m_w) + pt.m_Lf*H[j])*(H[i]*dens*detJ);
                vec3d kwJ = pt.m_aft*(-dens/(pt.m_ef+1)*H[i]*H[j]*detJ);
                
                ke[i7+3][j7  ] += Kwu(0,0); ke[i7+3][j7+1] += Kwu(0,1); ke[i7+3][j7+2] += Kwu(0,2);
                ke[i7+4][j7  ] += Kwu(1,0); ke[i7+4][j7+1] += Kwu(1,1); ke[i7+4][j7+2] += Kwu(1,2);
                ke[i7+5][j7  ] += Kwu(2,0); ke[i7+5][j7+1] += Kwu(2,1); ke[i7+5][j7+2] += Kwu(2,2);
                
                ke[i7+3][j7+3] += Kww(0,0); ke[i7+3][j7+4] += Kww(0,1); ke[i7+3][j7+5] += Kww(0,2);
                ke[i7+4][j7+3] += Kww(1,0); ke[i7+4][j7+4] += Kww(1,1); ke[i7+4][j7+5] += Kww(1,2);
                ke[i7+5][j7+3] += Kww(2,0); ke[i7+5][j7+4] += Kww(2,1); ke[i7+5][j7+5] += Kww(2,2);
                
                ke[i7+3][j7+6] += kwJ.x;
                ke[i7+4][j7+6] += kwJ.y;
                ke[i7+5][j7+6] += kwJ.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::Update(const FETimeInfo& tp)
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
void FEFluidFSIDomain3D::UpdateElementStress(int iel, const FETimeInfo& tp)
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
    
    // nodal coordinates
    const int NELN = FEElement::MAX_NODES;
    vec3d r0[NELN], r[NELN];
    vec3d vs[NELN];
    vec3d a[NELN];
    vec3d w[NELN];
    vec3d aw[NELN];
    double e[NELN];
    double ae[NELN];
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

		// FSI material point data
		ft.m_w = el.Evaluate(w, n);
		ft.m_Jdot = (Jt - Jp) / dt*dtrans;
		ft.m_aw = el.Evaluate(aw, n)*dtrans;

		// fluid material point data
		pt.m_efdot = el.Evaluate(ae, n)*dtrans;
		pt.m_vft = ept.m_v + ft.m_w;
		mat3d Gradw = Gradient(el, w, n);
		mat3d Lw = Gradw*Fi;
		pt.m_Lf = ept.m_L + Lw;
        pt.m_ef = el.Evaluate(e, n);
		vec3d Gradef = Gradient(el, e, n);
		pt.m_gradef = Fi.transpose()*Gradef;

		// fluid acceleration
		pt.m_aft = ept.m_a + ft.m_aw + pt.m_Lf*ft.m_w;

        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, tp);
        
		// calculate the stresses at this material point
		pt.m_sf = m_pMat->Fluid()->Stress(mp);
        ft.m_ss = m_pMat->Solid()->Stress(mp);
        ept.m_s = pt.m_sf + ft.m_ss;

		// calculate the fluid pressure
		pt.m_pf = m_pMat->Fluid()->Pressure(pt.m_ef);
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::InertialForces(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
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
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::ElementInertialForce(FESolidElement& el, vector<double>& fe)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, n;
    
    // jacobian determinant
    double detJ;
    
    const double* H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*	gw = el.GaussWeights();
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        double dens = m_pMat->Fluid()->Density(mp);
        
        // calculate the jacobian
        detJ = detJt(el, n, tp.alphaf)*gw[n];
        
        H = el.H(n);
        
        for (i=0; i<neln; ++i)
        {
            vec3d f = pt.m_aft*(dens*H[i]*detJ);
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[7*i+3] -= f.x;
            fe[7*i+4] -= f.y;
            fe[7*i+5] -= f.z;
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSIDomain3D::Serialize(DumpStream& ar)
{
	FESolidDomain::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_pMat;
	ar & m_sseps;
    ar & m_dofU & m_dofV & m_dofW & m_dofAW;
    ar & m_dofSU & m_dofR;
    ar & m_dof;
    ar & m_dofEF & m_dofAEF;
}
