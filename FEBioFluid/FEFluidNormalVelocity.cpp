//
//  FEFluidNormalVelocity.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEFluidNormalVelocity.h"
#include "FECore/FEModel.h"
#include "FECore/FEElemElemList.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/log.h"
#include "FECore/LinearSolver.h"

//=============================================================================
BEGIN_PARAMETER_LIST(FEFluidNormalVelocity, FESurfaceLoad)
ADD_PARAMETER(m_velocity, FE_PARAM_DOUBLE    , "velocity");
ADD_PARAMETER(m_VC      , FE_PARAM_DATA_ARRAY, "value"   );
ADD_PARAMETER(m_bpv     , FE_PARAM_BOOL      , "prescribe_nodal_velocities");
ADD_PARAMETER(m_bpar    , FE_PARAM_BOOL      , "parabolic"  );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEFluidNormalVelocity::FEFluidNormalVelocity(FEModel* pfem) : FESurfaceLoad(pfem), m_VC(FE_DOUBLE)
{
    m_velocity = 0.0;
    m_VC.set(1.0);
    m_bpv = true;
    m_bpar = false;
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    m_dofE = pfem->GetDOFIndex("e");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidNormalVelocity::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
    m_VC.Create(ps);
}

//-----------------------------------------------------------------------------
void FEFluidNormalVelocity::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int N = el.Nodes();
    lm.resize(N);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        lm[i] = id[m_dofE];
    }
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the prescribed normal velocity
void FEFluidNormalVelocity::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
    FEModel& fem = R.GetFEModel();
    
    vector<double> fe;
    vector<int> elm;
    
    vec3d r0[FEElement::MAX_NODES];
    
    int i, n;
    int npr = (int)m_VC.size();
    for (int iel=0; iel<npr; ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        int ndof = el.Nodes();
        fe.resize(ndof);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (i=0; i<neln; ++i) r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;
        
        double* Gr, *Gs;
        double* N;
        double* w  = el.GaussWeights();
        
        vec3d dxr, dxs;
        
        // get the velocity at the integration point
        double vn;
        
        // repeat over integration points
        zero(fe);
        for (n=0; n<nint; ++n)
        {
            N  = el.H(n);
            Gr = el.Gr(n);
            Gs = el.Gs(n);
            
            // calculate the tangent vectors
            vn = 0;
            dxr = dxs = vec3d(0,0,0);
            for (i=0; i<neln; ++i)
            {
                vn += N[i]*m_VN[el.m_lnode[i]];
                dxr += r0[i]*Gr[i];
                dxs += r0[i]*Gs[i];
            }
            
            vn *= m_velocity;
            double da = (dxr ^ dxs).norm();
            
            for (i=0; i<neln; ++i)
                fe[i] += N[i]*vn*w[n]*da;
        }
        
        // get the element's LM vector and adjust it
        UnpackLM(el, elm);
        
        // add element force vector to global force vector
        R.Assemble(el.m_node, elm, fe);
    }
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidNormalVelocity::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    FEMesh* mesh = ps->GetMesh();

    // evaluate surface normals
    vector<vec3d> sn(ps->Elements(),vec3d(0,0,0));
    m_nu.resize(ps->Nodes(),vec3d(0,0,0));
    m_VN.resize(ps->Nodes(),0);
    vector<int> nf(ps->Nodes(),0);
    vec3d r0[FEElement::MAX_NODES];
    for (int iel=0; iel<ps->Elements(); ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int i=0; i<neln; ++i) {
            r0[i] = mesh->Node(el.m_node[i]).m_r0;
            m_VN[el.m_lnode[i]] += m_VC.get<double>(iel);
            ++nf[el.m_lnode[i]];
        }
        
        double* Nr, *Ns;
        double* w  = el.GaussWeights();
        
        vec3d dxr, dxs;
        
        // repeat over integration points
        for (int n=0; n<nint; ++n)
        {
            Nr = el.Gr(n);
            Ns = el.Gs(n);
            
            // calculate the tangent vectors at integration point
            dxr = dxs = vec3d(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                dxr += r0[i]*Nr[i];
                dxs += r0[i]*Ns[i];
            }
            
            sn[iel] = dxr ^ dxs;
            sn[iel].unit();
        }
        
        // evaluate nodal normals by averaging surface normals
        for (int i=0; i<neln; ++i)
            m_nu[el.m_lnode[i]] += sn[iel];
    }
    
    for (int i=0; i<ps->Nodes(); ++i) {
        m_nu[i].unit();
        m_VN[i] /= nf[i];
    }
    
    // Set parabolic velocity profile if requested.
    // This will override velocity boundary cards in m_VC
    // and nodal cards in m_VN
    if (m_bpar) return SetParabolicVelocity();
    
    return true;
}

//-----------------------------------------------------------------------------
//! Mark the nodes with prescribed velocities
void FEFluidNormalVelocity::MarkVelocity()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();
    
    int id;
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        id = node.m_ID[m_dofVX]; if (id >= 0) node.m_ID[m_dofVX] = -id-2;
        id = node.m_ID[m_dofVY]; if (id >= 0) node.m_ID[m_dofVY] = -id-2;
        id = node.m_ID[m_dofVZ]; if (id >= 0) node.m_ID[m_dofVZ] = -id-2;
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the velocities
void FEFluidNormalVelocity::SetVelocity()
{
    // prescribe this velocity at the nodes
    FESurface* ps = &GetSurface();
    
    int id;
    for (int i=0; i<ps->Nodes(); ++i)
    {
        // evaluate the velocity
        vec3d v = m_nu[i]*(m_velocity*m_VN[i]);
        FENode& node = ps->Node(i);
        if (node.m_ID[m_dofVX] < -1) node.set(m_dofVX, v.x);
        if (node.m_ID[m_dofVY] < -1) node.set(m_dofVY, v.y);
        if (node.m_ID[m_dofVZ] < -1) node.set(m_dofVZ, v.z);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate normal velocities by solving Poiseuille flow across the surface
bool FEFluidNormalVelocity::SetParabolicVelocity()
{
    // find surface boundary nodes
    FEElemElemList EEL;
    FESurface* ps = &GetSurface();
    EEL.Create(ps);

    vector<bool> boundary(ps->Nodes(), false);
    for (int i=0; i<ps->Elements(); ++i) {
        FESurfaceElement& el = ps->Element(i);
        for (int j=0; j<el.facet_edges(); ++j) {
            FEElement* nel = EEL.Neighbor(i, j);
            if (nel == nullptr) {
                int en[3] = {-1,-1,-1};
                el.facet_edge(j, en);
                boundary[en[0]] = true;
                boundary[en[1]] = true;
                if (en[2] > -1) boundary[en[2]] = true;
            }
        }
    }
    
    // count number of non-boundary nodes
    int neq = 0;
    vector<int> glm(ps->Nodes(),-1);
    for (int i=0; i<ps->Nodes(); ++i)
        if (!boundary[i]) glm[i] = neq++;
    if (neq == 0)
    {
        felog.printbox("FATAL ERROR","Unable to set parabolic fluid normal velocity\n");
        return false;
    }
    
    // create a linear solver
    LinearSolver*		plinsolve;	//!< the linear solver
    FEGlobalMatrix*		pK;			//!< stiffness matrix
    FECoreKernel& fecore = FECoreKernel::GetInstance();
    plinsolve = fecore.CreateLinearSolver(SKYLINE_SOLVER);
    if (plinsolve == 0)
    {
        felog.printbox("FATAL ERROR","Unknown solver type selected\n");
        return false;
    }
    
    SparseMatrix* pS = plinsolve->CreateSparseMatrix(REAL_SYMMETRIC);
    pK = new FEGlobalMatrix(pS);
    if (pK == 0)
    {
        felog.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
        return false;
    }
    // build matrix profile for normal velocity at non-boundary nodes
    pK->build_begin(neq);
    for (int i=0; i<ps->Elements(); ++i) {
        FESurfaceElement& el = ps->Element(i);
        vector<int> elm(el.Nodes(),-1);
        for (int j=0; j<el.Nodes(); ++j)
            elm[j] = glm[el.m_lnode[j]];
        pK->build_add(elm);
    }
    pK->build_end();
    
    // create global vector
    vector<double> v;           //!< normal velocity solution
    vector<double> rhs;         //!< right-hand-side
    vector<double> Fr;          //!< reaction forces
    v.assign(neq, 0);
    rhs.assign(neq, 0);
    Fr.assign(neq, 0);
    FEModel pfem;
    FEGlobalVector pR(pfem, rhs, Fr);
    
    // calculate the global matrix and vector
    matrix ke;
    vector<double> fe;
    vector<int> lm;
    
    for (int m=0; m<ps->Elements(); ++m)
    {
        // get the surface element
        FESurfaceElement& el = ps->Element(m);
        
        int neln = el.Nodes();

        // get the element stiffness matrix
        ke.resize(neln, neln);
        lm.resize(neln);
        fe.resize(neln);
        vector<vec3d> gradN(neln);
        
        // calculate stiffness
        int nint = el.GaussPoints();
        
        // gauss weights
        double* w = el.GaussWeights();
        
        // nodal coordinates
        FEMesh& mesh = *ps->GetMesh();
        vec3d rt[FEElement::MAX_NODES];
        for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;
        
        // repeat over integration points
        ke.zero();
        zero(fe);
        vec3d gcnt[2];
        for (int n=0; n<nint; ++n)
        {
            double* N = el.H(n);
            double* Gr = el.Gr(n);
            double* Gs = el.Gs(n);
            ps->ContraBaseVectors(el, n, gcnt);
            
            vec3d dxr(0,0,0), dxs(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                dxr += rt[i]*Gr[i];
                dxs += rt[i]*Gs[i];
                gradN[i] = gcnt[0]*Gr[i] + gcnt[1]*Gs[i];
            }
            
            double da = (dxr ^ dxs).norm();
            
            // calculate stiffness component
            for (int i=0; i<neln; ++i) {
                fe[i] += N[i]*w[n]*da;
                for (int j=0; j<neln; ++j)
                    ke[i][j] += (gradN[i]*gradN[j])*w[n]*da;
            }
        }
        
        // get the element's LM vector
        for (int j=0; j<el.Nodes(); ++j)
            lm[j] = glm[el.m_lnode[j]];
        
        // assemble element matrix in global stiffness matrix
        pK->Assemble(ke, lm);
        pR.Assemble(lm, fe);
    }

    // solve linear system
    plinsolve->PreProcess();
    plinsolve->Factor();
    if (plinsolve->BackSolve(v, rhs) == false)
    {
        felog.printbox("FATAL ERROR","Unable to solve for parabolic fluid normal velocity\n");
        return false;
    }
    
    // set the nodal normal velocity scale factors
    for (int i=0; i<ps->Nodes(); ++i) {
        if (glm[i] == -1) m_VN[i] = 0;
        else m_VN[i] = v[glm[i]];
    }
    
    // evaluate net area and volumetric flow rate
    double A = 0, Q = 0;
    for (int m=0; m<ps->Elements(); ++m)
    {
        // get the surface element
        FESurfaceElement& el = ps->Element(m);
        
        int neln = el.Nodes();
        int nint = el.GaussPoints();
        double* w = el.GaussWeights();
        
        // nodal coordinates
        FEMesh& mesh = *ps->GetMesh();
        vec3d rt[FEElement::MAX_NODES];
        for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;
        
        // repeat over integration points
        for (int n=0; n<nint; ++n)
        {
            double* N = el.H(n);
            double* Gr = el.Gr(n);
            double* Gs = el.Gs(n);

            double vn = 0;
            vec3d dxr(0,0,0), dxs(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                vn += N[i]*m_VN[el.m_lnode[i]];
                dxr += rt[i]*Gr[i];
                dxs += rt[i]*Gs[i];
            }
            
            double da = (dxr ^ dxs).norm();
            
            for (int i=0; i<neln; ++i) {
                A += N[i]*w[n]*da;
                Q += N[i]*vn*w[n]*da;
            }
        }
    }
    
    // normalize nodal velocity cards
    double vbar = Q/A;
    for (int i=0; i<ps->Nodes(); ++i) m_VN[i] /= vbar;
    
    return true;
}
