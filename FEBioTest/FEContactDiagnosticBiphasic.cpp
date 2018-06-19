//
//  FEContactDiagnosticBiphasic.cpp
//  FEBioTest
//
//  Created by Gerard Ateshian on 5/1/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FEContactDiagnosticBiphasic.h"
#include "FEBioMix/FEBiphasicSolver.h"
#include "FEBioMix/FEBiphasicSolidDomain.h"
#include "FEBioMix/FESlidingInterfaceBiphasic.h"
#include "FEBioMech/FEResidualVector.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactDiagnosticBiphasic::FEContactDiagnosticBiphasic(FEModel& fem) : FEDiagnostic(fem)
{
    m_pscn = 0;
    
    FEAnalysis* pstep = new FEAnalysis(&fem);
    
    // create a new solver
    FESolver* pnew_solver = fecore_new<FESolver>(FESOLVER_ID, "biphasic", &fem);
    assert(pnew_solver);
    pnew_solver->m_bsymm = false;
    pstep->SetFESolver(pnew_solver);
    
    fem.AddStep(pstep);
    fem.SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEContactDiagnosticBiphasic::~FEContactDiagnosticBiphasic()
{
    
}

//-----------------------------------------------------------------------------
// Helper function to print a matrix
void FEContactDiagnosticBiphasic::print_matrix(matrix& m)
{
    int i, j;
    int N = m.rows();
    int M = m.columns();
    
    felog.printf("\n    ");
    for (i=0; i<N; ++i) felog.printf("%15d ", i);
    felog.printf("\n----");
    for (i=0; i<N; ++i) felog.printf("----------------", i);
    
    for (i=0; i<N; ++i)
    {
        felog.printf("\n%2d: ", i);
        for (j=0; j<M; ++j)
        {
            felog.printf("%15lg ", m[i][j]);
        }
    }
    felog.printf("\n");
}

//-----------------------------------------------------------------------------
// Helper function to print a sparse matrix
void FEContactDiagnosticBiphasic::print_matrix(SparseMatrix& m)
{
    int i, j;
	int N = m.Rows();
    int M = m.Columns();
    
    felog.printf("\n    ");
    for (i=0; i<N; ++i) felog.printf("%15d ", i);
    felog.printf("\n----");
    for (i=0; i<N; ++i) felog.printf("----------------", i);
    
    for (i=0; i<N; ++i)
    {
        felog.printf("\n%2d: ", i);
        for (j=0; j<M; ++j)
        {
            felog.printf("%15lg ", m.get(i,j));
        }
    }
    felog.printf("\n");
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEContactDiagnosticBiphasic::Init()
{
    if (m_pscn == 0) return false;
    
    if (m_pscn->Init() == false) return false;
    
    return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
bool FEContactDiagnosticBiphasic::Run()
{
    // get the solver
    FEModel& fem = GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
	fem.GetTime().timeIncrement = pstep->m_dt0 = dt;
    pstep->m_tstart = 0;
    pstep->m_tend = dt;
    pstep->m_final_time = dt;
    pstep->Activate();
    FEBiphasicSolver& solver = static_cast<FEBiphasicSolver&>(*pstep->GetFESolver());
    solver.m_bsymm = false;
    solver.Init();
    
    // make sure contact data is up to data
    solver.UpdateContact();
    
    // create the stiffness matrix
    solver.CreateStiffness(true);
    
    // get the stiffness matrix
    FEGlobalMatrix& K = solver.GetStiffnessMatrix();
    SparseMatrix& K0 = *K;
    
    // build the stiffness matrix
    K.Zero();
    solver.ContactStiffness();
    
    print_matrix(K0);
    
    // calculate the derivative of the residual
    matrix K1;
    deriv_residual(K1);
    
    print_matrix(K1);
    
    // calculate difference matrix
    int ndpn = 4;
    int N = mesh.Nodes();
    const int ndof = N*ndpn;
    matrix Kd; Kd.resize(ndof, ndof);
    double kij, kmax = 0, k0;
    int i0=-1, j0=-1;
    for (int i=0; i<ndof; ++i)
        for (int j=0; j<ndof; ++j)
        {
            Kd(i,j) = K1(i,j) - K0.get(i,j);
            k0 = fabs(K0.get(i,j));
            if (k0 < 1e-15) k0 = 1;
            kij = 100.0*fabs(Kd(i,j))/k0;
            if (kij > kmax)
            {
                kmax = kij;
                i0 = i;
                j0 = j;
            }
        }
    
    print_matrix(Kd);
    
    felog.printf("\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);
    
    return (kmax < 1e-3);
}

//-----------------------------------------------------------------------------
void FEContactDiagnosticBiphasic::deriv_residual(matrix& K)
{
    // get the solver
    FEModel& fem = GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
	fem.GetTime().timeIncrement = pstep->m_dt0 = dt;
    pstep->m_tstart = 0;
    pstep->m_tend = dt;
    pstep->m_final_time = dt;
    FEBiphasicSolver& solver = static_cast<FEBiphasicSolver&>(*pstep->GetFESolver());
    
    // get the DOFs
    const int dof_x = fem.GetDOFIndex("x");
    const int dof_y = fem.GetDOFIndex("y");
    const int dof_z = fem.GetDOFIndex("z");
    const int dof_p = fem.GetDOFIndex("p");
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    int ndpn = 4;
    int N = mesh.Nodes();
    int ndof = N*ndpn;
    
    solver.UpdateContact();
    
    // first calculate the initial residual
    vector<double> R0; R0.assign(ndof, 0);
    vector<double> dummy(R0);
    FEResidualVector RHS0(fem, R0, dummy);
    solver.ContactForces(RHS0);
    
    // now calculate the perturbed residuals
    K.resize(ndof,ndof);
    int i, j, nj;
    double dx = 1e-8;
    vector<double> R1(ndof);
    for (j=0; j<ndof; ++j)
    {
        FENode& node = mesh.Node(j/ndpn);
        nj = j%ndpn;
        
        switch (nj)
        {
            case 0: node.inc(dof_x, dx); node.m_rt.x += dx; break;
            case 1: node.inc(dof_y, dx); node.m_rt.y += dx; break;
            case 2: node.inc(dof_z, dx); node.m_rt.z += dx; break;
            case 3: node.inc(dof_p, dx); break;
        }
        
        solver.UpdateContact();
        
        zero(R1);
        FEResidualVector RHS1(fem, R1, dummy);
        solver.ContactForces(RHS1);
        
        switch (nj)
        {
            case 0: node.dec(dof_x, dx); node.m_rt.x -= dx; break;
            case 1: node.dec(dof_y, dx); node.m_rt.y -= dx; break;
            case 2: node.dec(dof_z, dx); node.m_rt.z -= dx; break;
            case 3: node.dec(dof_p, dx); break;
        }
        
        solver.UpdateContact();
        
        for (i=0; i<ndof; ++i) K(i,j) = (R0[i] - R1[i])/dx;
    }
}

//-----------------------------------------------------------------------------
FEDiagnosticScenario* FEContactDiagnosticBiphasic::CreateScenario(const std::string& sname)
{
    if (sname == "hex8") return (m_pscn = new FEContactBiphasicTangentHex8(this));
    else if (sname == "hex20") return (m_pscn = new FEContactBiphasicTangentHex20(this));
    return 0;
}

//////////////////////////////////////////////////////////////////////
// Biphasic Contact Tangent Diagnostic for hex8 Elements
//////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEContactBiphasicTangentHex8, FEContactBiphasicScenario)
ADD_PARAMETER(m_dt      , FE_PARAM_DOUBLE, "time_step"     );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEContactBiphasicTangentHex8::FEContactBiphasicTangentHex8(FEDiagnostic* pdia) : FEContactBiphasicScenario(pdia)
{
    m_dt = 1;
}

//-----------------------------------------------------------------------------
bool FEContactBiphasicTangentHex8::Init()
{
    vec3d r[16] = {
        vec3d(0,0,0), vec3d(1,0,0), vec3d(1,1,0), vec3d(0,1,0),
        vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1),
        vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1),
        vec3d(0,0,2), vec3d(1,0,2), vec3d(1,1,2), vec3d(0,1,2)
    };
    
    double p[16] = {
        1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1
    };
    
    FEModel& fem = GetDiagnostic()->GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // --- create the geometry ---
    int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
    
    // get the DOFs
    const int dof_p = fem.GetDOFIndex("p");
    
    // currently we simply assume a two-element contact problem
    // so we create two elements
    const double eps = 0.5;
    mesh.CreateNodes(16);
    mesh.SetDOFS(MAX_DOFS);
    vec3d dr(0,0,eps);
    for (int i=0; i<8; ++i) {
        mesh.Node(i).m_r0 = r[i];
        mesh.Node(i).set(dof_p, p[i]);
    }
    for (int i=8; i<16; ++i) {
        mesh.Node(i).m_r0 = r[i] - dr;
        mesh.Node(i).set(dof_p, p[i]+eps);
    }
    
    for (int i=0; i<16; ++i)
    {
        FENode& node = mesh.Node(i);
        node.m_rt = node.m_r0;
    }
    
    // get the material
    FEMaterial* pm = fem.GetMaterial(0);
    
    // get the one-and-only domain
    FEBiphasicSolidDomain* pbd = new FEBiphasicSolidDomain(&fem);
    pbd->SetMaterial(pm);
    pbd->Create(2, FE_HEX8G8);
	pbd->SetMatID(0);
    mesh.AddDomain(pbd);
    
    FESolidElement& el0 = pbd->Element(0);
    FESolidElement& el1 = pbd->Element(1);
    
    el0.SetID(1);
    el0.m_node[0] = 0;
    el0.m_node[1] = 1;
    el0.m_node[2] = 2;
    el0.m_node[3] = 3;
    el0.m_node[4] = 4;
    el0.m_node[5] = 5;
    el0.m_node[6] = 6;
    el0.m_node[7] = 7;
    
    el1.SetID(2);
    el1.m_node[0] = 8;
    el1.m_node[1] = 9;
    el1.m_node[2] = 10;
    el1.m_node[3] = 11;
    el1.m_node[4] = 12;
    el1.m_node[5] = 13;
    el1.m_node[6] = 14;
    el1.m_node[7] = 15;
    
    pbd->CreateMaterialPointData();
    
    // --- create the sliding interface ---
    FESlidingInterfaceBiphasic* ps = new FESlidingInterfaceBiphasic(&fem);
    ps->m_atol = 0.1;
    ps->m_epsn = 1;
    ps->m_epsp = 1;
    ps->m_btwo_pass = false;
    ps->m_nsegup = 0;
    ps->m_bautopen = true;
    ps->m_bsymm = false;
    ps->m_stol = 0.01;
    FESlidingSurfaceBiphasic& ms = ps->m_ms;
    ms.Create(1, FE_QUAD4G4);
    ms.Element(0).m_node[0] = 4;
    ms.Element(0).m_node[1] = 5;
    ms.Element(0).m_node[2] = 6;
    ms.Element(0).m_node[3] = 7;
    FESlidingSurfaceBiphasic& ss = ps->m_ss;
    ss.Create(1, FE_QUAD4G4);
    ss.Element(0).m_node[0] = 11;
    ss.Element(0).m_node[1] = 10;
    ss.Element(0).m_node[2] = 9;
    ss.Element(0).m_node[3] = 8;
    fem.AddSurfacePairConstraint(ps);
    
    return true;
}

//////////////////////////////////////////////////////////////////////
// Biphasic Contact Tangent Diagnostic for hex20 Elements
//////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEContactBiphasicTangentHex20, FEContactBiphasicScenario)
ADD_PARAMETER(m_dt      , FE_PARAM_DOUBLE, "time_step"     );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEContactBiphasicTangentHex20::FEContactBiphasicTangentHex20(FEDiagnostic* pdia) : FEContactBiphasicScenario(pdia)
{
    m_dt = 1;
}

//-----------------------------------------------------------------------------
bool FEContactBiphasicTangentHex20::Init()
{
    vec3d r[40] = {
        vec3d(0,0,0), vec3d(0,0,0.5), vec3d(0,0,1), vec3d(0,0.5,0),
        vec3d(0,0.5,1), vec3d(0,1,0), vec3d(0,1,0.5), vec3d(0,1,1),
        vec3d(0.5,0,0), vec3d(0.5,0,1), vec3d(0.5,1,0), vec3d(0.5,1,1),
        vec3d(1,0,0), vec3d(1,0,0.5), vec3d(1,0,1), vec3d(1,0.5,0),
        vec3d(1,0.5,1), vec3d(1,1,0), vec3d(1,1,0.5), vec3d(1,1,1),
        vec3d(0,0,1), vec3d(0,0,1.5), vec3d(0,0,2), vec3d(0,0.5,1),
        vec3d(0,0.5,2), vec3d(0,1,1), vec3d(0,1,1.5), vec3d(0,1,2),
        vec3d(0.5,0,1), vec3d(0.5,0,2), vec3d(0.5,1,1), vec3d(0.5,1,2),
        vec3d(1,0,1), vec3d(1,0,1.5), vec3d(1,0,2), vec3d(1,0.5,1),
        vec3d(1,0.5,2), vec3d(1,1,1), vec3d(1,1,1.5), vec3d(1,1,2),
    };
    
    double p[40] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    };
    
    int el0n[20] = {
        1,    13,    18,     6,     3,    15,    20,     8,     9,    16,
        11,     4,    10,    17,    12,     5,     2,    14,    19,     7
    };
    
    int el1n[20] = {
        21,    33,    38,    26,    23,    35,    40,    28,    29,    36,
        31,    24,    30,    37,    32,    25,    22,    34,    39,    27
    };
    
    int msn[8] = {
        21,    26,    38,    33,    24,    31,    36,    29
    };
    
    int ssn[9] = {
        3,    15,    20,     8,    10,    17,    12,     5
    };
    
    FEModel& fem = GetDiagnostic()->GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // --- create the geometry ---
    int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
    
    // get the DOFs
    const int dof_p = fem.GetDOFIndex("p");
    
    // currently we simply assume a two-element contact problem
    // so we create two elements
    const double eps = 0.5;
    mesh.CreateNodes(40);
    mesh.SetDOFS(MAX_DOFS);
    vec3d dr(0,0,eps);
    for (int i=0; i<20; ++i) {
        mesh.Node(i).m_r0 = r[i];
        mesh.Node(i).set(dof_p, p[i]);
    }
    for (int i=20; i<40; ++i) {
        mesh.Node(i).m_r0 = r[i] - dr;
        mesh.Node(i).set(dof_p, p[i]+eps);
    }
    
    for (int i=0; i<40; ++i)
    {
        FENode& node = mesh.Node(i);
        node.m_rt = node.m_r0;
    }
    
    // get the material
    FEMaterial* pm = fem.GetMaterial(0);
    
    // get the one-and-only domain
    FEBiphasicSolidDomain* pbd = new FEBiphasicSolidDomain(&fem);
    pbd->SetMaterial(pm);
    pbd->Create(2, FE_HEX20G27);
	pbd->SetMatID(0);
    mesh.AddDomain(pbd);
    
    FESolidElement& el0 = pbd->Element(0);
    FESolidElement& el1 = pbd->Element(1);
    
    el0.SetID(1);
    for (int i=0; i<20; ++i) {
        el0.m_node[i] = el0n[i] - 1;
    }
    
    el1.SetID(2);
    for (int i=0; i<20; ++i) {
        el1.m_node[i] = el1n[i] - 1;
    }
    
    pbd->CreateMaterialPointData();
    
    // --- create the sliding interface ---
    FESlidingInterfaceBiphasic* ps = new FESlidingInterfaceBiphasic(&fem);
    ps->m_atol = 0.1;
    ps->m_epsn = 1;
    ps->m_epsp = 1;
    ps->m_btwo_pass = false;
    ps->m_nsegup = 0;
    ps->m_bautopen = true;
    ps->m_bsymm = false;
    FESlidingSurfaceBiphasic& ms = ps->m_ms;
    ms.Create(1, FE_QUAD8G9);
    for (int i=0; i<8; ++i)
        ms.Element(0).m_node[i] = msn[i] - 1;
    FESlidingSurfaceBiphasic& ss = ps->m_ss;
    ss.Create(1, FE_QUAD8G9);
    for (int i=0; i<8; ++i)
        ss.Element(0).m_node[i] = ssn[i] - 1;
    fem.AddSurfacePairConstraint(ps);
    
    return true;
}
