//
//  FEFluidDomain2D.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 12/15/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEFluidDomain2D.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "NumCore/LUSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEFluidDomain2D::FEFluidDomain2D(FEModel* pfem) : FEDomain2D(&pfem->GetMesh()), FEFluidDomain(pfem)
{
    m_pMat = 0;
    m_btrans = true;

    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
    m_dofEF  = pfem->GetDOFIndex("ef");
    
    m_dofWXP = pfem->GetDOFIndex("wxp");
    m_dofWYP = pfem->GetDOFIndex("wyp");
    m_dofWZP = pfem->GetDOFIndex("wzp");
    m_dofEFP  = pfem->GetDOFIndex("efp");
    
    m_dofAWX = pfem->GetDOFIndex("awx");
    m_dofAWY = pfem->GetDOFIndex("awy");
    m_dofAWZ = pfem->GetDOFIndex("awz");
    m_dofAEF = pfem->GetDOFIndex("aef");
    
    m_dofAWXP = pfem->GetDOFIndex("awxp");
    m_dofAWYP = pfem->GetDOFIndex("awyp");
    m_dofAWZP = pfem->GetDOFIndex("awzp");
    m_dofAEFP = pfem->GetDOFIndex("aefp");

	// list the degrees of freedom
	// (This allows the FEBomain base class to handle several tasks such as UnpackLM)
	vector<int> dof;
	dof.push_back(m_dofWX);
	dof.push_back(m_dofWY);
	dof.push_back(m_dofEF);
	SetDOFList(dof);
}

//-----------------------------------------------------------------------------
// \todo I don't think this is being used
FEFluidDomain2D& FEFluidDomain2D::operator = (FEFluidDomain2D& d)
{
    m_Elem = d.m_Elem;
    m_pMesh = d.m_pMesh;
    return (*this);
}

//-----------------------------------------------------------------------------
//! Assign material
void FEFluidDomain2D::SetMaterial(FEMaterial* pmat)
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
//! \todo The material point initialization needs to move to the base class.
bool FEFluidDomain2D::Init()
{
    // initialize base class
	FEDomain2D::Init();
    
    // get the elements material
    FEFluid* pme = m_pMat;
    
    // assign local coordinate system to each integration point
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FEElement2D& el = m_Elem[i];
        for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
    }
    
    // check for initially inverted elements
    int ninverted = 0;
    for (int i=0; i<Elements(); ++i)
    {
        FEElement2D& el = Element(i);
        
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
void FEFluidDomain2D::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    const int NE = FEElement::MAX_NODES;
    vec3d x0[NE], r0, v;
    FEMesh& m = *GetMesh();
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FEElement2D& el = m_Elem[i];
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
                felog.printbox("ERROR", "Negative jacobian was detected.");
                throw DoRunningRestart();
            }
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain2D::InternalForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEElement2D& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 3*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        //#pragma omp critical
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEFluidDomain2D::ElementInternalForce(FEElement2D& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[2][2], detJ;
    
    mat3ds sv;
    vec3d gradp;
    
    const double *H, *Gr, *Gs;
    
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
        
        vec3d g1(Ji[0][0],Ji[0][1],0.0);
        vec3d g2(Ji[1][0],Ji[1][1],0.0);
        
        // get the viscous stress tensor for this integration point
        sv = m_pMat->GetViscous()->Stress(mp);
        // get the gradient of the elastic pressure
        gradp = pt.m_gradJf*m_pMat->Tangent_Pressure_Strain(mp);
        
        H = el.H(n);
        Gr = el.Hr(n);
        Gs = el.Hs(n);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i];
        
        // Jdot/J
        double dJoJ = pt.m_Jfdot/pt.m_Jf;
        
        for (i=0; i<neln; ++i)
        {
            vec3d fs = sv*gradN[i] + gradp*H[i];
            double fJ = dJoJ*H[i] + gradN[i]*pt.m_vft;
            
            // calculate internal force
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[3*i  ] -= fs.x*detJ;
            fe[3*i+1] -= fs.y*detJ;
            fe[3*i+2] -= fJ*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain2D::BodyForce(FEGlobalVector& R, const FETimeInfo& tp, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEElement2D& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 3*el.Nodes();
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

void FEFluidDomain2D::ElementBodyForce(FEBodyForce& BF, FEElement2D& el, vector<double>& fe)
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
            fe[3*i  ] -= H[i]*dens*f.x*detJ;
            fe[3*i+1] -= H[i]*dens*f.y*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEFluidDomain2D::ElementBodyForceStiffness(FEBodyForce& BF, FEElement2D &el, matrix &ke)
{
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
        detJ = detJ0(el, n)*gw[n];
        
        H = el.H(n);
        
        double dens = m_pMat->Density(mp);
        
        // get the force
        f = BF.force(mp);
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i) {
            for (int j=0; j<neln; ++j)
            {
                k = f*(-H[i]*H[j]*dens/pt.m_Jf*detJ);
                ke[ndof*i  ][ndof*j+2] += k.x;
                ke[ndof*i+1][ndof*j+2] += k.y;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEFluidDomain2D::ElementMaterialStiffness(FEElement2D &el, matrix &ke)
{
    int i, i3, j, j3, n;

	double dt = GetFEModel()->GetTime().timeIncrement;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double *H, *Gr, *Gs;
    
    // jacobian
    double Ji[2][2], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    // get materials
    FEViscousFluid* m_pViscous = m_pMat->GetViscous();
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],0.0);
        vec3d g2(Ji[1][0],Ji[1][1],0.0);
        
        H = el.H(n);
        Gr = el.Hr(n);
        Gs = el.Hs(n);
        
        // setup the material point
        // NOTE: deformation gradient and determinant have already been evaluated in the stress routine
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
        // get the tangents
        double dpdJ = m_pMat->Tangent_Pressure_Strain(mp);
        mat3ds svJ = m_pViscous->Tangent_Strain(mp);
        tens4ds cv = m_pViscous->Tangent_RateOfDeformation(mp);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i];
        
        // evaluate stiffness matrix
        for (i=0, i3=0; i<neln; ++i, i3 += 3)
        {
            for (j=0, j3 = 0; j<neln; ++j, j3 += 3)
            {
                mat3d Kv = vdotTdotv(gradN[i], cv, gradN[j])*detJ;
                vec3d kv = (pt.m_gradJf*H[j] - gradN[j]*pt.m_Jf)*(H[i]*detJ);
                vec3d kJ = (mat3dd(-dpdJ) + svJ)*gradN[i]*(H[j]*detJ);
                double k = (H[j]*((1.0*m_btrans)/dt - pt.m_Lf.trace()) + gradN[j]*pt.m_vft)*(H[i]*detJ);
                
                ke[i3  ][j3  ] += Kv(0,0);
                ke[i3  ][j3+1] += Kv(0,1);
                ke[i3  ][j3+2] += kJ.x;
                
                ke[i3+1][j3  ] += Kv(1,0);
                ke[i3+1][j3+1] += Kv(1,1);
                ke[i3+1][j3+2] += kJ.y;
                
                ke[i3+2][j3  ] += kv.x;
                ke[i3+2][j3+1] += kv.y;
                ke[i3+2][j3+2] += k;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain2D::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FEElement2D& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 3*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate material stiffness
        ElementMaterialStiffness(el, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain2D::MassMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FEElement2D& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 3*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementMassMatrix(el, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain2D::BodyForceStiffness(FESolver* psolver, const FETimeInfo& tp, FEBodyForce& bf)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(GetMaterial()); assert(pme);
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FEElement2D& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 3*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementBodyForceStiffness(bf, el, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function

void FEFluidDomain2D::ElementStiffness(int iel, matrix& ke)
{
    FEElement2D& el = Element(iel);
    
    // calculate material stiffness (i.e. constitutive component)
    ElementMaterialStiffness(el, ke);
    
}

//-----------------------------------------------------------------------------
//! calculates element inertial stiffness matrix
void FEFluidDomain2D::ElementMassMatrix(FEElement2D& el, matrix& ke)
{
    int i, i3, j, j3, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    // gradient of shape functions
    vector<vec3d> gradN(neln);
    
    double *H;
    double *Gr, *Gs;
    
    // jacobian
    double Ji[2][2], detJ;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        detJ = invjac0(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],0.0);
        vec3d g2(Ji[1][0],Ji[1][1],0.0);
        
        H = el.H(n);
        
        Gr = el.Hr(n);
        Gs = el.Hs(n);
        
        // setup the material point
        // NOTE: deformation gradient and determinant have already been evaluated in the stress routine
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
        
		double dt = GetFEModel()->GetTime().timeIncrement;
        double dens = m_pMat->Density(mp);
        
        // evaluate spatial gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i];
        
        // evaluate stiffness matrix
        for (i=0, i3=0; i<neln; ++i, i3 += 3)
        {
            for (j=0, j3 = 0; j<neln; ++j, j3 += 3)
            {
                mat3d Mv = ((mat3dd(1)*(m_btrans/dt) + pt.m_Lf)*H[j] + mat3dd(gradN[j]*pt.m_vft))*(H[i]*dens*detJ);
                vec3d mJ = pt.m_aft*(-H[i]*H[j]*dens/pt.m_Jf*detJ);
                
                ke[i3  ][j3  ] += Mv(0,0);
                ke[i3  ][j3+1] += Mv(0,1);
                ke[i3  ][j3+2] += mJ.x;
                
                ke[i3+1][j3  ] += Mv(1,0);
                ke[i3+1][j3+1] += Mv(1,1);
                ke[i3+1][j3+2] += mJ.y;
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain2D::Update(const FETimeInfo& tp)
{
    // TODO: This is temporary hack for running micro-materials in parallel.
    //	     Evaluating the stress for a micro-material will make FEBio solve
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
                if (NegativeJacobian::m_boutput) e.print();
            }
        }
    }
    
    // reset the logfile mode
    felog.SetMode(nmode);
    
    // if we encountered an error, we request a running restart
    if (berr)
    {
        if (NegativeJacobian::m_boutput == false) felog.printbox("ERROR", "Negative jacobian was detected.");
        throw DoRunningRestart();
    }
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
void FEFluidDomain2D::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double alphaf = tp.alphaf;
    double alpham = tp.alpham;
    
    // get the solid element
    FEElement2D& el = m_Elem[iel];
    
    // get the number of integration points
    int nint = el.GaussPoints();
    
    // number of nodes
    int neln = el.Nodes();
    
    // nodal coordinates
    vec3d vt[FEElement::MAX_NODES], vp[FEElement::MAX_NODES];
    vec3d at[FEElement::MAX_NODES], ap[FEElement::MAX_NODES];
    double et[FEElement::MAX_NODES], ep[FEElement::MAX_NODES];
    double aet[FEElement::MAX_NODES], aep[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) {
        FENode& node = m_pMesh->Node(el.m_node[j]);
        vt[j] = node.get_vec3d(m_dofWX, m_dofWY, m_dofWZ);
        vp[j] = node.get_vec3d(m_dofWXP, m_dofWYP, m_dofWZP);
        at[j] = node.get_vec3d(m_dofAWX, m_dofAWY, m_dofAWZ);
        ap[j] = node.get_vec3d(m_dofAWXP, m_dofAWYP, m_dofAWZP);
        et[j] = node.get(m_dofEF);
        ep[j] = node.get(m_dofEFP);
        aet[j] = node.get(m_dofAEF);
        aep[j] = node.get(m_dofAEFP);
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
        pt.m_Jf = 1 + el.Evaluate(et, n)*alphaf + el.Evaluate(ep, n)*(1-alphaf);
        pt.m_gradJf = gradient(el, et, n)*alphaf + gradient(el, ep, n)*(1-alphaf);
        pt.m_Jfdot = pt.m_gradJf*pt.m_vft;
        if (m_btrans) pt.m_Jfdot += el.Evaluate(aet, n)*alpham + el.Evaluate(aep, n)*(1-alpham);
        
        // calculate the stress at this material point
        pt.m_sf = m_pMat->Stress(mp);
        
        // calculate the fluid pressure
        pt.m_pf = m_pMat->Pressure(mp);
    }
}

//-----------------------------------------------------------------------------
void FEFluidDomain2D::InertialForces(FEGlobalVector& R, const FETimeInfo& tp)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEElement2D& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 3*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInertialForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        //#pragma omp critical
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEFluidDomain2D::ElementInertialForce(FEElement2D& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian determinant
    double detJ;
    
    mat3ds s;
    
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
            fe[3*i  ] -= f.x*detJ;
            fe[3*i+1] -= f.y*detJ;
        }
    }
}
