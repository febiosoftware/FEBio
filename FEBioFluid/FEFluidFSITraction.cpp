//
//  FEFluidFSITraction.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/14/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEFluidFSITraction.h"
#include "FECore/FEModel.h"
#include "FEFluid.h"
#include "FEFluidFSI.h"

//-----------------------------------------------------------------------------
// Parameter block for traction loads
BEGIN_PARAMETER_LIST(FEFluidFSITraction, FESurfaceLoad)
ADD_PARAMETER(m_s    , FE_PARAM_DOUBLE, "scale");
ADD_PARAMETER(m_bself, FE_PARAM_BOOL  , "self" );
END_PARAMETER_LIST()

//-----------------------------------------------------------------------------
//! constructor
FEFluidFSITraction::FEFluidFSITraction(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_K = 0;
    m_s = 1.0;
    m_bself = false;

    // get the degrees of freedom
    m_dofX = pfem->GetDOFIndex("x");
    m_dofY = pfem->GetDOFIndex("y");
    m_dofZ = pfem->GetDOFIndex("z");
    m_dofWX = pfem->GetDOFIndex("wx");
    m_dofWY = pfem->GetDOFIndex("wy");
    m_dofWZ = pfem->GetDOFIndex("wz");
    m_dofEF = pfem->GetDOFIndex("ef");
    m_dofEFP = pfem->GetDOFIndex("efp");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidFSITraction::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidFSITraction::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    if (!m_bself) ps->SetInterfaceStatus(true);
    ps->Init();
    
    // get the list of fluid-FSI elements connected to this interface
    FEModel* fem = GetFEModel();
    FEMesh* mesh = ps->GetMesh();
    int NF = ps->Elements();
    m_elem.resize(NF);
    for (int j=0; j<NF; ++j)
    {
        FESurfaceElement& el = ps->Element(j);
        // extract the first of two elements on this interface
        m_elem[j] = mesh->FindElementFromID(el.m_elem[0]);
        // get its material and check if FluidFSI
        FEMaterial* pm = fem->GetMaterial(m_elem[j]->GetMatID());
        FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(pm);
        if (pfsi == nullptr) {
            // extract the second of two elements on this interface
            m_elem[j] = mesh->FindElementFromID(el.m_elem[1]);
            pm = fem->GetMaterial(m_elem[j]->GetMatID());
            pfsi = dynamic_cast<FEFluidFSI*>(pm);
            if (pfsi == nullptr) return false;
        }
    }
    
    // get the fluid bulk modulus to calculate pressure from dilatation
    // get material data from first surface element
    // assuming the entire surface interfaces between two same domains
    FESurfaceElement& el = ps->Element(0);
    if (!m_bself) {
        if ((el.m_elem[0] == -1) || (el.m_elem[1] == -1)) return false;
        FEMesh* mesh = ps->GetMesh();
        FEElement* pe1 = mesh->FindElementFromID(el.m_elem[0]);
        FEElement* pe2 = mesh->FindElementFromID(el.m_elem[1]);
        if ((pe1 == nullptr) || (pe2 == nullptr)) return false;
        // get the material
        FEMaterial* pm1 = GetFEModel()->GetMaterial(pe1->GetMatID());
        FEMaterial* pm2 = GetFEModel()->GetMaterial(pe2->GetMatID());
        FEFluidFSI* fsi = dynamic_cast<FEFluidFSI*> (pm1);
        if (fsi == nullptr) fsi = dynamic_cast<FEFluidFSI*> (pm2);
        if (fsi == nullptr) return false;
        // get the bulk modulus
        m_K = fsi->Fluid()->m_k;
    }
    else {
        if (el.m_elem[0] == -1) return false;
        FEMesh* mesh = ps->GetMesh();
        FEElement* pe1 = mesh->FindElementFromID(el.m_elem[0]);
        if (pe1 == nullptr) return false;
        // get the material
        FEMaterial* pm1 = GetFEModel()->GetMaterial(pe1->GetMatID());
        FEFluidFSI* fsi = dynamic_cast<FEFluidFSI*> (pm1);
        if (fsi == nullptr) return false;
        // get the bulk modulus
        m_K = fsi->Fluid()->m_k;
    }

    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = *GetSurface().GetMesh();
    int N = el.Nodes();
    lm.resize(N*7);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        lm[7*i  ] = id[m_dofX];
        lm[7*i+1] = id[m_dofY];
        lm[7*i+2] = id[m_dofZ];
        lm[7*i+3] = id[m_dofWX];
        lm[7*i+4] = id[m_dofWY];
        lm[7*i+5] = id[m_dofWZ];
        lm[7*i+6] = id[m_dofEF];
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
    FESurface& surf = GetSurface();
    int npr = surf.Elements();
#pragma omp parallel for
    for (int i=0; i<npr; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        FESurfaceElement& el = m_psurf->Element(i);
        
        // calculate nodal normal tractions
        int neln = el.Nodes();
        vector<double> tn(neln);
        
        int ndof = 7*neln;
        fe.resize(ndof);
        
        ElementForce(el, fe, tp, i);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // add element force vector to global force vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the element force

void FEFluidFSITraction::ElementForce(FESurfaceElement& el, vector<double>& fe, const FETimeInfo& tp, const int iel)
{
    // nr integration points
    int nint = el.GaussPoints();
    
    // nr of element nodes
    int neln = el.Nodes();
    
    FEModel* fem = GetFEModel();
    
    // nodal coordinates
    FEMesh& mesh = *m_psurf->GetMesh();
    vec3d rt[FEElement::MAX_NODES];
    double et[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) {
        FENode& node = mesh.Node(el.m_node[j]);
        rt[j] = node.m_rt*tp.alpha + node.m_rp*(1-tp.alpha);
        et[j] = node.get(m_dofEF)*tp.alphaf + node.get(m_dofEFP)*(1-tp.alphaf);
    }
    
    // Get the fluid stress from the fluid-FSI element
    mat3ds sv(mat3dd(0));
    FEElement* pe = m_elem[iel];
    int pint = pe->GaussPoints();
    FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(fem->GetMaterial(pe->GetMatID()));
    for (int n=0; n<pint; ++n)
    {
        FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
        sv += pfsi->Fluid()->GetViscous()->Stress(mp);
    }
    sv /= pint;
    
    // repeat over integration points
    zero(fe);
    double* w  = el.GaussWeights();
    for (int j=0; j<nint; ++j)
    {
        double* N  = el.H(j);
        double* Gr = el.Gr(j);
        double* Gs = el.Gs(j);
        
        // evaluate fluid dilatation and covariant basis vectors
        // at integration point
        double ef = 0;
        vec3d gr(0,0,0), gs(0,0,0);
        for (int i=0; i<neln; ++i) {
            ef += et[i]*N[i];
            gr += rt[i]*Gr[i];
            gs += rt[i]*Gs[i];
        }
        
        vec3d gt = gr ^ gs;
        
        vec3d f = (sv*gt + gt*(m_K*ef))*(-m_s*w[j]);
        
        for (int i=0; i<neln; ++i)
        {
            fe[7*i  ] += N[i]*f.x;
            fe[7*i+1] += N[i]*f.y;
            fe[7*i+2] += N[i]*f.z;
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
{
    FESurface& surf = GetSurface();
    int npr = surf.Elements();
#pragma omp parallel for
    for (int m=0; m<npr; ++m)
    {
        matrix ke;
        vector<int> lm;
        
        // get the surface element
        FESurfaceElement& el = m_psurf->Element(m);
        
        // calculate nodal normal tractions
        int neln = el.Nodes();
        vector<double> tn(neln);
        
        // get the element stiffness matrix
        int ndof = 7*neln;
        ke.resize(ndof, ndof);
        
        // calculate pressure stiffness
        ElementStiffness(el, ke, tp, m);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEFluidFSITraction::ElementStiffness(FESurfaceElement& el, matrix& ke, const FETimeInfo& tp, const int iel)
{
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    mat3dd I(1);
    
    FEModel* fem = GetFEModel();
    FESurface* ps = &GetSurface();

    // gauss weights
    double* w = el.GaussWeights();
    
    double dt = tp.timeIncrement;
    double alpha = tp.alpha;
    double a = tp.gamma/(tp.beta*dt);

    vector<vec3d> gradN(neln);

    // nodal coordinates
    FEMesh& mesh = *m_psurf->GetMesh();
    vec3d rt[FEElement::MAX_NODES];
    double et[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) {
        FENode& node = mesh.Node(el.m_node[j]);
        rt[j] = node.m_rt*tp.alpha + node.m_rp*(1-tp.alpha);
        et[j] = node.get(m_dofEF)*tp.alphaf + node.get(m_dofEFP)*(1-tp.alphaf);
    }
    
    // Get the fluid stress and its tangents from the fluid-FSI element
    mat3ds sv(mat3dd(0)), svJ(mat3dd(0));
    tens4ds cv; cv.zero();
    mat3d Ls; Ls.zero();
    FEElement* pe = m_elem[iel];
    int pint = pe->GaussPoints();
    FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(fem->GetMaterial(pe->GetMatID()));
    for (int n=0; n<pint; ++n)
    {
        FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
        FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
        sv += pfsi->Fluid()->GetViscous()->Stress(mp);
        svJ += pfsi->Fluid()->GetViscous()->Tangent_Strain(mp);
        cv += pfsi->Fluid()->Tangent_RateOfDeformation(mp);
        Ls += ep.m_L;
    }
    sv /= pint;
    svJ /= pint;
    cv /= pint;
    Ls /= pint;
    mat3d M = mat3dd(a) - Ls;
    
    // repeat over integration points of surface element
    ke.zero();
    for (int k=0; k<nint; ++k)
    {
        double* N = el.H(k);
        double* Gr = el.Gr(k);
        double* Gs = el.Gs(k);
        
        // evaluate fluid dilatation and covariant basis vectors
        double ef = 0;
        vec3d gr(0,0,0), gs(0,0,0);
        for (int i=0; i<neln; ++i) {
            ef += et[i]*N[i];
            gr += rt[i]*Gr[i];
            gs += rt[i]*Gs[i];
        }

        // evaluate fluid pressure
        double p = m_K*ef*w[k]*m_s;
        
        vec3d gt = gr ^ gs;
        vec3d f = gt*(-m_K*w[k]*m_s);

        vec3d gcnt[2], gcntp[2];
        ps->ContraBaseVectors(el, k, gcnt);
        ps->ContraBaseVectorsP(el, k, gcntp);
        for (int i=0; i<neln; ++i)
            gradN[i] = (gcnt[0]*alpha + gcntp[0]*(1-alpha))*Gr[i] +
            (gcnt[1]*alpha + gcntp[1]*(1-alpha))*Gs[i];
        
        // calculate stiffness component
        int i, j, i7, j7;
        for (i=0, i7=0; i<neln; ++i, i7+=7)
            for (j=0, j7=0; j<neln; ++j, j7+=7)
            {
                vec3d v = gr*Gs[j] - gs*Gr[j];
                mat3d A; A.skew(v);
                mat3d Kv = vdotTdotv(gt, cv, gradN[j]);
                
                mat3d Kuu = (sv*A + Kv*M)*(-N[i]*w[k]*m_s) - A*(N[i]*p); Kuu *= alpha;
                mat3d Kuw = Kv*(-N[i]*w[k]*m_s); Kuw *= alpha;
                vec3d kuJ = svJ*gt*(-N[i]*N[j]*w[k]*m_s) + f*(N[i]*N[j]); kuJ *= alpha;
                
                ke[i7  ][j7  ] -= Kuu(0,0); ke[i7  ][j7+1] -= Kuu(0,1); ke[i7  ][j7+2] -= Kuu(0,2);
                ke[i7+1][j7  ] -= Kuu(1,0); ke[i7+1][j7+1] -= Kuu(1,1); ke[i7+1][j7+2] -= Kuu(1,2);
                ke[i7+2][j7  ] -= Kuu(2,0); ke[i7+2][j7+1] -= Kuu(2,1); ke[i7+2][j7+2] -= Kuu(2,2);
                
                ke[i7  ][j7+3] -= Kuw(0,0); ke[i7  ][j7+4] -= Kuw(0,1); ke[i7  ][j7+5] -= Kuw(0,2);
                ke[i7+1][j7+3] -= Kuw(1,0); ke[i7+1][j7+4] -= Kuw(1,1); ke[i7+1][j7+5] -= Kuw(1,2);
                ke[i7+2][j7+3] -= Kuw(2,0); ke[i7+2][j7+4] -= Kuw(2,1); ke[i7+2][j7+5] -= Kuw(2,2);
                
                ke[i7  ][j7+6] -= kuJ.x;
                ke[i7+1][j7+6] -= kuJ.y;
                ke[i7+2][j7+6] -= kuJ.z;
            }
    }
}

//-----------------------------------------------------------------------------

void FEFluidFSITraction::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
}
