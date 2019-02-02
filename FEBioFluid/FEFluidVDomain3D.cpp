//
//  FEFluidVDomain3D.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 1/30/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEFluidVDomain3D.h"
#include "FEFluidVSolver.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEFluidVDomain3D::FEFluidVDomain3D(FEModel* pfem) : FESolidDomain(pfem), FEFluidDomain(pfem)
{
    m_pMat = 0;
    m_btrans = true;
    
    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
    
    m_dofWXP = pfem->GetDOFIndex("wxp");
    m_dofWYP = pfem->GetDOFIndex("wyp");
    m_dofWZP = pfem->GetDOFIndex("wzp");
    
    m_dofAWX = pfem->GetDOFIndex("awx");
    m_dofAWY = pfem->GetDOFIndex("awy");
    m_dofAWZ = pfem->GetDOFIndex("awz");
    
    m_dofAWXP = pfem->GetDOFIndex("awxp");
    m_dofAWYP = pfem->GetDOFIndex("awyp");
    m_dofAWZP = pfem->GetDOFIndex("awzp");
    
    // list the degrees of freedom
    // (This allows the FEBomain base class to handle several tasks such as UnpackLM)
    vector<int> dof;
    dof.push_back(m_dofWX);
    dof.push_back(m_dofWY);
    dof.push_back(m_dofWZ);
    SetDOFList(dof);
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEFluidVDomain3D& FEFluidVDomain3D::operator = (FEFluidVDomain3D& d)
{
    m_Elem = d.m_Elem;
    m_pMesh = d.m_pMesh;
    return (*this);
}

//-----------------------------------------------------------------------------
//! Assign material
void FEFluidVDomain3D::SetMaterial(FEMaterial* pmat)
{
    FEDomain::SetMaterial(pmat);
    if (pmat)
    {
        m_pMat = dynamic_cast<FEFluidV*>(pmat);
        assert(m_pMat);
    }
    else m_pMat = 0;
}

//-----------------------------------------------------------------------------
//! \todo The material point initialization needs to move to the base class.
bool FEFluidVDomain3D::Init()
{
    // initialize base class
    FESolidDomain::Init();
    
    // check for initially inverted elements
    int ninverted = 0;
    for (int i=0; i<Elements(); ++i)
    {
        FESolidElement& el = Element(i);
        
        int nint = el.GaussPoints();
        for (int n=0; n<nint; ++n)
        {
            double J0 = detJ0(el, n);
            if (J0 <= 0)
            {
                felog.printf("**************************** E R R O R ****************************\n");
                felog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.GetID());
                felog.printf("Jacobian = %lg\n", J0);
                felog.printf("Did you use the right node numbering?\n");
                felog.printf("Nodes:");
                for (int l=0; l<el.Nodes(); ++l)
                {
                    felog.printf("%d", el.m_node[l]+1);
                    if (l+1 != el.Nodes()) felog.printf(","); else felog.printf("\n");
                }
                felog.printf("*******************************************************************\n\n");
                ++ninverted;
            }
        }
    }
    
    return (ninverted == 0);
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEFluidVDomain3D::PreSolveUpdate(const FETimeInfo& timeInfo)
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
            FEFluidVMaterialPoint& vpt = *mp.ExtractData<FEFluidVMaterialPoint>();
            pt.m_r0 = el.Evaluate(x0, j);
            vpt.m_Jfp = vpt.m_Jft;
            vpt.m_dJfp = vpt.m_dJft;

            if (pt.m_Jf <= 0) {
                felog.printbox("ERROR", "Negative jacobian was detected.");
                throw DoRunningRestart();
            }
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidVDomain3D::InternalForces(FEGlobalVector& R, const FETimeInfo& tp)
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
        int ndof = 3*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForce(el, fe, tp);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        //#pragma omp critical
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEFluidVDomain3D::ElementInternalForce(FESolidElement& el, vector<double>& fe, const FETimeInfo& tp)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJ;
    
    mat3ds s, dsJ;
    
    const double *H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double*    gw = el.GaussWeights();
    
    double fJ = 0, kJJ = 0;
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);

    // evaluate fJ and kJJ by integrating over the entire element
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n];
        fJ += (pt.m_Lf.trace() - pt.m_Jfdot/pt.m_Jf)*detJ;
        kJJ += (pt.m_Jfdot/pt.m_Jf - ksi/dt)/pt.m_Jf*detJ;
    }

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
        s = pt.m_sf;
        dsJ = m_pMat->Fluid()->Tangent_Strain(mp);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        for (i=0; i<neln; ++i)
        {
            vec3d kvJ = dsJ*gradN[i];
            vec3d fv = s*gradN[i];
            vec3d f = fv;
            if (kJJ != 0) f -= kvJ*(fJ/kJJ);
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[3*i  ] -= f.x*detJ;
            fe[3*i+1] -= f.y*detJ;
            fe[3*i+2] -= f.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidVDomain3D::BodyForce(FEGlobalVector& R, const FETimeInfo& tp, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 3*el.Nodes();
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

void FEFluidVDomain3D::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe, const FETimeInfo& tp)
{
    double *H;
    double* gw = el.GaussWeights();
    vec3d b;
    
    // number of nodes and integration points
    int neln = el.Nodes();
    int nint = el.GaussPoints();

    // nodal coordinates
    vec3d r0[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
        r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
    
    // jacobian
    double detJ;
    double fJ = 0, kJJ = 0;
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);

    // evaluate fJ and kJJ by integrating over the entire element
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n];
        fJ += (pt.m_Lf.trace() - pt.m_Jfdot/pt.m_Jf)*detJ;
        kJJ += (pt.m_Jfdot/pt.m_Jf - ksi/dt)/pt.m_Jf*detJ;
    }
    
    // loop over integration points
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();
        double dens = m_pMat->Fluid()->Density(mp);
        
        pt.m_r0 = el.Evaluate(r0, n);
        
        detJ = detJ0(el, n)*gw[n];
        H = el.H(n);

        // get the force
        b = BF.force(mp);
        
        for (int i=0; i<neln; ++i)
        {
            vec3d fb = b*(H[i]*dens);
            vec3d kb = -fb/pt.m_Jf;
            vec3d f = fb;
            if (kJJ != 0) f -= kb*(fJ/kJJ);
            fe[3*i  ] -= f.x*detJ;
            fe[3*i+1] -= f.y*detJ;
            fe[3*i+2] -= f.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidVDomain3D::BodyForceStiffness(FESolver* psolver, const FETimeInfo& tp, FEBodyForce& bf)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(GetMaterial()); assert(pme);
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FESolidElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 3*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementBodyForceStiffness(bf, el, ke, tp);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEFluidVDomain3D::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke, const FETimeInfo& tp)
{
    int neln = el.Nodes();
    int nint = el.GaussPoints();
    int ndof = ke.columns()/neln;
    
    double *H, *Gr, *Gs, *Gt;
    double* gw = el.GaussWeights();
    vec3d b;
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    // jacobian
    double Ji[3][3], detJ;
    double kJJ = 0;
    vector<vec3d> kb(neln,vec3d(0,0,0));
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);
    
    // evaluate kJJ by integrating over the entire element
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        H = el.H(n);
        
        double dens = m_pMat->Fluid()->Density(mp);
        
        // get the force
        b = BF.force(mp);
        
        kJJ += (pt.m_Jfdot/pt.m_Jf - ksi/dt)/pt.m_Jf*detJ;
        for (int i=0; i<neln; ++i)
            kb[i] -= b*dens*H[i]/pt.m_Jf*detJ;
    }
    
    // loop over integration points
    for (int n=0; n<nint; ++n)
    {
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n];
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
        
        for (int i=0; i<neln; ++i) {
            for (int j=0; j<neln; ++j)
            {
                mat3d K(0,0,0,0,0,0,0,0,0);
                if (kJJ!= 0) K = (kb[i]/kJJ) & (gradN[j]*detJ);
                ke[ndof*i  ][ndof*j  ] -= K(0,0); ke[ndof*i  ][ndof*j+1] -= K(0,1); ke[ndof*i  ][ndof*j+2] -= K(0,2);
                ke[ndof*i+1][ndof*j  ] -= K(1,0); ke[ndof*i+1][ndof*j+1] -= K(1,1); ke[ndof*i+1][ndof*j+2] -= K(1,2);
                ke[ndof*i+2][ndof*j  ] -= K(2,0); ke[ndof*i+2][ndof*j+1] -= K(2,1); ke[ndof*i+2][ndof*j+2] -= K(2,2);
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidVDomain3D::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FESolidElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 3*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate material stiffness
        ElementStiffness(el, ke, tp);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEFluidVDomain3D::ElementStiffness(FESolidElement &el, matrix &ke, const FETimeInfo& tp)
{
    int i, i3, j, j3, n;
    
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
    
    // evaluate kJJ by integrating over the entire element
    double kJJ = 0;
    vector<vec3d> kvJ(neln,vec3d(0,0,0));
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        mat3ds dsJ = m_pMat->Fluid()->Tangent_Strain(mp);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i) {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
            kvJ[i] += dsJ*gradN[i]*detJ;
        }
        
        kJJ += (pt.m_Jfdot/pt.m_Jf - ksi/dt)/pt.m_Jf*detJ;
    }
    
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
        
        // get the tangents
        tens4ds cv = m_pMat->Fluid()->Tangent_RateOfDeformation(mp);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // evaluate stiffness matrix
        for (i=0, i3=0; i<neln; ++i, i3 += 3)
        {
            for (j=0, j3 = 0; j<neln; ++j, j3 += 3)
            {
                mat3d Kvv = vdotTdotv(gradN[i], cv, gradN[j]);
                mat3d K = Kvv;
                if (kJJ !=0) K -= (kvJ[i]/kJJ) & gradN[j];
                K *= detJ;
                ke[i3  ][j3  ] += K(0,0); ke[i3  ][j3+1] += K(0,1); ke[i3  ][j3+2] += K(0,2);
                ke[i3+1][j3  ] += K(1,0); ke[i3+1][j3+1] += K(1,1); ke[i3+1][j3+2] += K(1,2);
                ke[i3+2][j3  ] += K(2,0); ke[i3+2][j3+1] += K(2,1); ke[i3+2][j3+2] += K(2,2);
            }
        }
    }
}


//-----------------------------------------------------------------------------
void FEFluidVDomain3D::InertialForces(FEGlobalVector& R, const FETimeInfo& tp)
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
        int ndof = 3*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInertialForce(el, fe, tp);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        //#pragma omp critical
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
void FEFluidVDomain3D::ElementInertialForce(FESolidElement& el, vector<double>& fe, const FETimeInfo& tp)
{
    int i, n;
    
    // jacobian determinant
    double Ji[3][3], detJ;
    
    const double* H, *Gr, *Gs, *Gt;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double*    gw = el.GaussWeights();
    
    double fJ = 0, kJJ = 0;
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);
    
    // evaluate fJ and kJJ by integrating over the entire element
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // calculate the jacobian
        detJ = detJ0(el, n)*gw[n];
        fJ += (pt.m_Lf.trace() - pt.m_Jfdot/pt.m_Jf)*detJ;
        kJJ += (pt.m_Jfdot/pt.m_Jf - ksi/dt)/pt.m_Jf*detJ;
    }
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        double dens = m_pMat->Fluid()->Density(mp);
        
        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];

        for (i=0; i<neln; ++i)
        {
            vec3d fa = pt.m_aft*(H[i]*dens);
            vec3d mvJ = -fa/pt.m_Jf;
            vec3d f = fa;
            if (kJJ != 0) f -= mvJ*(fJ/kJJ);
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[3*i  ] -= f.x*detJ;
            fe[3*i+1] -= f.y*detJ;
            fe[3*i+2] -= f.z*detJ;
        }
    }
}
//-----------------------------------------------------------------------------
void FEFluidVDomain3D::MassMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FESolidElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 3*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementMassMatrix(el, ke, tp);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! calculates element inertial stiffness matrix
void FEFluidVDomain3D::ElementMassMatrix(FESolidElement& el, matrix& ke, const FETimeInfo& tp)
{
    int i, i3, j, j3, n;
    
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
    
    // evaluate kJJ by integrating over the entire element
    double kJJ = 0;
    vector<vec3d> mvJ(neln,vec3d(0,0,0));
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        double dens = m_pMat->Fluid()->Density(mp);
        H = el.H(n);
        kJJ += (pt.m_Jfdot/pt.m_Jf - ksi/dt)/pt.m_Jf*detJ;
        for (int i=0; i<neln; ++i)
            mvJ[i] -= pt.m_aft*(dens*H[i]/pt.m_Jf)*detJ;
    }
    
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
        for (i=0, i3=0; i<neln; ++i, i3 += 3)
        {
            for (j=0, j3 = 0; j<neln; ++j, j3 += 3)
            {
                mat3d Mvv = ((mat3dd(ksi*m_btrans/dt) + pt.m_Lf)*H[j] + mat3dd(gradN[j]*pt.m_vft))*(H[i]*dens*detJ);
                mat3d K = Mvv;
                if (kJJ != 0) Mvv -= (mvJ[i]/kJJ) & gradN[j];
                ke[i3  ][j3  ] += K(0,0); ke[i3  ][j3+1] += K(0,1); ke[i3  ][j3+2] += K(0,2);
                ke[i3+1][j3  ] += K(1,0); ke[i3+1][j3+1] += K(1,1); ke[i3+1][j3+2] += K(1,2);
                ke[i3+2][j3  ] += K(2,0); ke[i3+2][j3+1] += K(2,1); ke[i3+2][j3+2] += K(2,2);
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidVDomain3D::Update(const FETimeInfo& tp)
{
    // TODO: This is temporary hack for running micro-materials in parallel.
    //         Evaluating the stress for a micro-material will make FEBio solve
    //       a new FE problem. We don't want to see the output of that problem.
    //       The logfile is a shared resource between the master FEM and the RVE
    //       in order not to corrupt the logfile we don't print anything for
    //       the RVE problem.
    // TODO: Maybe I need to create a new domain class for micro-material.
    Logfile::MODE nmode = felog.GetMode();
    felog.SetMode(Logfile::LOG_NEVER);
    
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
                felog.SetMode(nmode);
                berr = true;
                if (NegativeJacobian::DoOutput()) e.print();
            }
        }
    }
    
    // reset the logfile mode
    felog.SetMode(nmode);
    
    // if we encountered an error, we request a running restart
    if (berr)
    {
        if (NegativeJacobian::DoOutput() == false) felog.printbox("ERROR", "Negative jacobian was detected.");
        throw DoRunningRestart();
    }
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
void FEFluidVDomain3D::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double alphaf = tp.alphaf;
    double alpham = tp.alpham;
    double gamma = tp.gamma;
    
    // get the solid element
    FESolidElement& el = m_Elem[iel];
    
    // get the number of integration points
    int nint = el.GaussPoints();
    
    // number of nodes
    int neln = el.Nodes();

    // gradient of shape functions
    vector<vec3d> gradN(neln);

    double* gw = el.GaussWeights();

    // jacobian
    double Ji[3][3], detJ;
    const double *Gr, *Gs, *Gt;
    double fJ = 0, kJJ = 0;
    vector<vec3d> kJv(neln,vec3d(0,0,0));
    double dt = tp.timeIncrement;
    double ksi = tp.alpham/(tp.gamma*tp.alphaf);
    
    // evaluate fJ and kJJ by integrating over the entire element
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // calculate the jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        fJ += (pt.m_Lf.trace() - pt.m_Jfdot/pt.m_Jf)*detJ;
        kJJ += (pt.m_Jfdot/pt.m_Jf - ksi/dt)/pt.m_Jf*detJ;
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        // evaluate spatial gradient of shape functions
        for (int i=0; i<neln; ++i) {
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
            kJv[i] += gradN[i]*detJ;
        }
    }
    
    // nodal coordinates
    const int NELN = FEElement::MAX_NODES;
    vec3d vt[NELN], vp[NELN];
    vec3d at[NELN], ap[NELN];
    for (int j=0; j<neln; ++j) {
        FENode& node = m_pMesh->Node(el.m_node[j]);
        vt[j] = node.get_vec3d(m_dofWX, m_dofWY, m_dofWZ);
        vp[j] = node.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP);
        at[j] = node.get_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ);
        ap[j] = node.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
    }
    
    // evaluate the increment in dilatation
    double dJ = 0;
    if (kJJ != 0) {
        dJ = -fJ/kJJ;
        for (int j=0; j<neln; ++j) {
            dJ -= (kJv[j]*(vt[j] - vp[j]))/kJJ;
        }
    }
    
    // loop over the integration points and update
    // velocity, velocity gradient, acceleration
    // stress and pressure at the integration point
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        FEFluidVMaterialPoint& vpt = *(mp.ExtractData<FEFluidVMaterialPoint>());

        // material point data
        pt.m_vft = el.Evaluate(vt, n)*alphaf + el.Evaluate(vp, n)*(1-alphaf);
        pt.m_Lf = gradient(el, vt, n)*alphaf + gradient(el, vp, n)*(1-alphaf);
        pt.m_aft = pt.m_Lf*pt.m_vft;
        if (m_btrans) pt.m_aft += el.Evaluate(at, n)*alpham + el.Evaluate(ap, n)*(1-alpham);
        vpt.m_Jft = vpt.m_Jfp + dJ;
        vpt.m_dJft = vpt.m_dJfp*(1.-1./gamma) + (vpt.m_Jft - vpt.m_Jfp)/(gamma*dt);
        pt.m_Jf = alphaf*vpt.m_Jft + (1-alphaf)*vpt.m_Jfp;
        pt.m_Jfdot = 0;
        if (m_btrans) pt.m_Jfdot = pt.m_Jf*pt.m_Lf.trace();
//        if (m_btrans) pt.m_Jfdot = alpham*vpt.m_dJft + (1-alpham)*vpt.m_dJft;

        // calculate the stress at this material point
        pt.m_sf = m_pMat->Fluid()->Stress(mp);
        
        // calculate the fluid pressure
        pt.m_pf = m_pMat->Fluid()->Pressure(mp);
    }
}
