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
#include "FEFluidDomain3D.h"
#include "FEFluidSolver.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioFluid.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEFluidDomain3D::FEFluidDomain3D(FEModel* pfem) : FESolidDomain(pfem), FEFluidDomain(pfem), m_dofW(pfem), m_dofAW(pfem), m_dof(pfem)
{
    m_pMat = 0;
    m_btrans = true;

	if (pfem)
	{
		m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
		m_dofAW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_ACCELERATION));

		m_dofEF = pfem->GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
		m_dofAEF = pfem->GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION_TDERIV), 0);

		FEDofList dofs(pfem);
		dofs.AddDofs(m_dofW);
		dofs.AddDof(m_dofEF);
		m_dof = dofs;
	}
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEFluidDomain3D& FEFluidDomain3D::operator = (FEFluidDomain3D& d)
{
    m_Elem = d.m_Elem;
    m_pMesh = d.m_pMesh;
    return (*this);
}

//-----------------------------------------------------------------------------
//! serialize data to archive
void FEFluidDomain3D::Serialize(DumpStream& ar)
{
    FESolidDomain::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofW & m_dofAW & m_dof;
    ar & m_dofEF & m_dofAEF;
    ar & m_pMat;
}

//-----------------------------------------------------------------------------
// get total dof list
const FEDofList& FEFluidDomain3D::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
//! Assign material
void FEFluidDomain3D::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
    if (pmat)
    {
        m_pMat = dynamic_cast<FEFluid*>(pmat);
        assert(m_pMat);
    }
    else m_pMat = 0;
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEFluidDomain3D::PreSolveUpdate(const FETimeInfo& timeInfo)
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
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain3D::InternalForces(FEGlobalVector& R)
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
        int ndof = 4*el.Nodes();
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

void FEFluidDomain3D::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJ;
    
    mat3ds sv;
    vec3d gradp;
    
    const double *H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double*	gw = el.GaussWeights();

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
        sv = m_pMat->GetViscous()->Stress(mp);
        // get the gradient of the elastic pressure
        gradp = pt.m_gradef*m_pMat->Tangent_Pressure_Strain(mp);
        
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
            vec3d fs = sv*gradN[i] + gradp*H[i];
            double fJ = dJoJ*H[i] + gradN[i]*pt.m_vft;
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[4*i  ] -= fs.x*detJ;
            fe[4*i+1] -= fs.y*detJ;
            fe[4*i+2] -= fs.z*detJ;
            fe[4*i+3] -= fJ*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain3D::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 4*el.Nodes();
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

void FEFluidDomain3D::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
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
            fe[4*i  ] -= H[i]*dens*f.x*detJ;
            fe[4*i+1] -= H[i]*dens*f.y*detJ;
            fe[4*i+2] -= H[i]*dens*f.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEFluidDomain3D::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int neln = el.Nodes();
    int ndof = ke.columns()/neln;
    
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
                ke[ndof*i  ][ndof*j+3] += k.x;
                ke[ndof*i+1][ndof*j+3] += k.y;
                ke[ndof*i+2][ndof*j+3] += k.z;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEFluidDomain3D::ElementStiffness(FESolidElement &el, matrix &ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i4, j, j4, n;
    
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
        mat3ds svJ = m_pMat->GetViscous()->Tangent_Strain(mp);
        tens4ds cv = m_pMat->Tangent_RateOfDeformation(mp);
        double dp = m_pMat->Tangent_Pressure_Strain(mp);
        double d2p = m_pMat->Tangent_Pressure_Strain_Strain(mp);
        // Jdot/J
        double dJoJ = pt.m_efdot/Jf;
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i4=0; i<neln; ++i, i4 += 4)
        {
            for (j=0, j4 = 0; j<neln; ++j, j4 += 4)
            {
                mat3d Kvv = vdotTdotv(gradN[i], cv, gradN[j]);
                vec3d kJv = (pt.m_gradef*(H[i]/Jf) + gradN[i])*H[j];
                vec3d kvJ = (svJ*gradN[i])*H[j] + (gradN[j]*dp+pt.m_gradef*(H[j]*d2p))*H[i];
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
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain3D::StiffnessMatrix(FELinearSystem& LS)
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
        int ndof = 4*el.Nodes();
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
void FEFluidDomain3D::MassMatrix(FELinearSystem& LS)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();

#pragma omp parallel for shared(NE)
    for (int iel=0; iel<NE; ++iel)
    {
		FESolidElement& el = m_Elem[iel];

        // element stiffness matrix
		FEElementMatrix ke(el);
        
        // create the element's stiffness matrix
        int ndof = 4*el.Nodes();
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
void FEFluidDomain3D::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
    for (int iel=0; iel<NE; ++iel)
    {
		FESolidElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);
        
        // create the element's stiffness matrix
        int ndof = 4*el.Nodes();
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
//! calculates element inertial stiffness matrix
void FEFluidDomain3D::ElementMassMatrix(FESolidElement& el, matrix& ke)
{
    const FETimeInfo& tp = GetFEModel()->GetTime();
    int i, i4, j, j4, n;
    
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
        
        double dens = m_pMat->Density(mp);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i4=0; i<neln; ++i, i4 += 4)
        {
            for (j=0, j4 = 0; j<neln; ++j, j4 += 4)
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
void FEFluidDomain3D::Update(const FETimeInfo& tp)
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
void FEFluidDomain3D::UpdateElementStress(int iel, const FETimeInfo& tp)
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
    double et[NELN], ep[NELN];
    double aet[NELN], aep[NELN];
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
    }
    
    // loop over the integration points and update
    // velocity, velocity gradient, acceleration
    // stress and pressure at the integration point
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // material point data
        pt.m_vft = el.Evaluate(vt, n)*alphaf + el.Evaluate(vp, n)*(1-alphaf);
        pt.m_Lf = gradient(el, vt, n)*alphaf + gradient(el, vp, n)*(1-alphaf);
        pt.m_aft = pt.m_Lf*pt.m_vft;
        if (m_btrans) pt.m_aft += el.Evaluate(at, n)*alpham + el.Evaluate(ap, n)*(1-alpham);
        pt.m_ef = el.Evaluate(et, n)*alphaf + el.Evaluate(ep, n)*(1-alphaf);
        pt.m_gradef = gradient(el, et, n)*alphaf + gradient(el, ep, n)*(1-alphaf);
        pt.m_efdot = pt.m_gradef*pt.m_vft;
        if (m_btrans) pt.m_efdot += el.Evaluate(aet, n)*alpham + el.Evaluate(aep, n)*(1-alpham);

        // calculate the stress at this material point
        pt.m_sf = m_pMat->Stress(mp);
        
        // calculate the fluid pressure
        pt.m_pf = m_pMat->Pressure(mp);
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain3D::InertialForces(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for shared(NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 4*el.Nodes();
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
void FEFluidDomain3D::ElementInertialForce(FESolidElement& el, vector<double>& fe)
{
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
        double dens = m_pMat->Density(mp);
        
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n];
        
        H = el.H(n);
        
        for (i=0; i<neln; ++i)
        {
            vec3d f = pt.m_aft*(dens*H[i]);
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[4*i  ] -= f.x*detJ;
            fe[4*i+1] -= f.y*detJ;
            fe[4*i+2] -= f.z*detJ;
        }
    }
}
