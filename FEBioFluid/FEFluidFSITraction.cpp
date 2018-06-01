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
//! constructor
FEFluidFSITraction::FEFluidFSITraction(FEModel* pfem) : FESurfaceLoad(pfem)
{
    // get the degrees of freedom
    m_dofX = pfem->GetDOFIndex("x");
    m_dofY = pfem->GetDOFIndex("y");
    m_dofZ = pfem->GetDOFIndex("z");
    m_dofSX = pfem->GetDOFIndex("sx");
    m_dofSY = pfem->GetDOFIndex("sy");
    m_dofSZ = pfem->GetDOFIndex("sz");
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
    ps->SetInterfaceStatus(true);
    ps->Init();
    
    // get the list of fluid-FSI elements connected to this interface
    FEModel* fem = GetFEModel();
    FEMesh* mesh = ps->GetMesh();
    int NF = ps->Elements();
    m_elem.resize(NF);
    m_K.resize(NF,0);
    m_s.resize(NF,1);
    m_bself.resize(NF, false);
    for (int j=0; j<NF; ++j)
    {
        FESurfaceElement& el = ps->Element(j);
        // extract the first of two elements on this interface
        m_elem[j] = mesh->FindElementFromID(el.m_elem[0]);
        if (el.m_elem[1] == -1) m_bself[j] = true;
        // get its material and check if FluidFSI
        FEMaterial* pm = fem->GetMaterial(m_elem[j]->GetMatID());
        FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(pm);
        if (pfsi) {
            m_K[j] = pfsi->Fluid()->m_k;
        }
        else if (!m_bself[j]) {
            // extract the second of two elements on this interface
            m_elem[j] = mesh->FindElementFromID(el.m_elem[1]);
            pm = fem->GetMaterial(m_elem[j]->GetMatID());
            pfsi = dynamic_cast<FEFluidFSI*>(pm);
            if (pfsi == nullptr) return false;
            m_s[j] = -1;
            m_K[j] = pfsi->Fluid()->m_k;
        }
        else
            return false;
    }

    // TODO: Deal with the case when the surface is a shell domain separating two FSI domains
    // that use different fluid bulk moduli
    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = *GetSurface().GetMesh();
    FESurfaceElement& fel = dynamic_cast<FESurfaceElement&>(el);
    FEElement* pe = mesh.FindElementFromID(fel.m_elem[0]);
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    FEFluidFSI* fsi = dynamic_cast<FEFluidFSI*> (pm);
    if (fsi == nullptr) pe = mesh.FindElementFromID(fel.m_elem[1]);
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

    // substitute interface dofs for solid-shell interfaces
    FESolidElement& sel = static_cast<FESolidElement&>(*pe);
    for (int i = 0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = mesh.Node(sel.m_node[i]);
            vector<int>& id = node.m_ID;
            int j = el.FindNode(node.GetID()-1);
            
            // first the displacement dofs
            lm[7*j  ] = id[m_dofSX];
            lm[7*j+1] = id[m_dofSY];
            lm[7*j+2] = id[m_dofSZ];
        }
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
        
        vec3d f = (sv*gt + gt*(m_K[iel]*ef))*(-m_s[iel]*w[j]);
        
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
        double p = m_K[iel]*ef*w[k]*m_s[iel];
        
        vec3d gt = gr ^ gs;
        vec3d f = gt*(-m_K[iel]*w[k]*m_s[iel]);

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
                
                mat3d Kuu = (sv*A + Kv*M)*(-N[i]*w[k]*m_s[iel]) - A*(N[i]*p); Kuu *= alpha;
                mat3d Kuw = Kv*(-N[i]*w[k]*m_s[iel]); Kuw *= alpha;
                vec3d kuJ = svJ*gt*(-N[i]*N[j]*w[k]*m_s[iel]) + f*(N[i]*N[j]); kuJ *= alpha;
                
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
