//
//  FEFluidTangentDiagnostic.cpp
//  FEBioTest
//
//  Created by Gerard Ateshian on 11/17/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#include "FEFluidTangentDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEBioFluid/FEFluidSolver.h"
#include "FEBioFluid/FEFluidDomain3D.h"
#include "FECore/log.h"
#include <FECore/FEPrescribedDOF.h>
#include <FECore/LoadCurve.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidTangentUniaxial, FEFluidScenario)
	ADD_PARAMETER(m_velocity, "fluid_velocity");
	ADD_PARAMETER(m_dt      , "time_step"     );
END_FECORE_CLASS();

BEGIN_FECORE_CLASS(FEFluidTangentUniaxialSS, FEFluidScenario)
	ADD_PARAMETER(m_velocity, "fluid_velocity");
	ADD_PARAMETER(m_dt      , "time_step"     );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFluidTangentUniaxial::FEFluidTangentUniaxial(FEDiagnostic* pdia) : FEFluidScenario(pdia)
{
    m_velocity = 0;
    m_dt = 1;
}

//-----------------------------------------------------------------------------
// Build the uniaxial loading scenario
// Cube with uniaxial velocity prescribed along x on left face and dilatation
// fixed on right face. Dynamic analysis.
bool FEFluidTangentUniaxial::Init()
{
    int i;
    vec3d r[8] = {
        vec3d(0,0,0), vec3d(1,0,0), vec3d(1,1,0), vec3d(0,1,0),
        vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1)
    };
    
    int BC[8][4] = {
        { 0,-1,-1, 0},{ 0,-1,-1,-1},{ 0,-1,-1,-1}, { 0,-1,-1, 0},
        { 0,-1,-1, 0},{ 0,-1,-1,-1},{ 0,-1,-1,-1}, { 0,-1,-1, 0}
    };
    
    FEModel& fem = GetDiagnostic()->GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    pstep->m_nanalysis = FE_DYNAMIC;

    int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	const int dof_WX = fem.GetDOFIndex("wx");
	const int dof_WY = fem.GetDOFIndex("wy");
	const int dof_WZ = fem.GetDOFIndex("wz");
	const int dof_EF  = fem.GetDOFIndex("ef");

    // --- create the FE problem ---
    // create the mesh
    FEMesh& m = fem.GetMesh();
    m.CreateNodes(8);
	m.SetDOFS(MAX_DOFS);
    for (i=0; i<8; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_r0 = r[i];
        n.m_rid = -1;
        
        // set displacement BC's
        if (BC[i][0] == -1) fem.AddFixedBC(i, dof_WX);
        if (BC[i][1] == -1) fem.AddFixedBC(i, dof_WY);
        if (BC[i][2] == -1) fem.AddFixedBC(i, dof_WZ);
        if (BC[i][3] == -1) fem.AddFixedBC(i, dof_EF);
    }
    
    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    
    // create a fluid domain
    FEFluidDomain3D* pd = new FEFluidDomain3D(&fem);
    pd->SetMaterial(pmat);
    pd->Create(1, FE_HEX8G8);
	pd->SetMatID(0);
    m.AddDomain(pd);
    FESolidElement& el = pd->Element(0);
    el.SetID(1);
    for (i=0; i<8; ++i) el.m_node[i] = i;
    
    pd->CreateMaterialPointData();
    
    // Add a loadcurve
	FELoadCurve* plc = new FELoadCurve(new FELinearFunction(&fem, 1.0, 0.0));
    fem.AddLoadCurve(plc);
    
    // Add a prescribed BC
    int nd[4] = {0, 3, 4, 7};
	FEPrescribedDOF* pdc = new FEPrescribedDOF(&fem);
    fem.AddPrescribedBC(pdc);
    pdc->SetDOF(dof_WX).SetScale(m_velocity, 0);
    for (i = 0; i<4; ++i) pdc->AddNode(nd[i]);
    
    return true;
}

//-----------------------------------------------------------------------------
FEFluidTangentUniaxialSS::FEFluidTangentUniaxialSS(FEDiagnostic* pdia) : FEFluidScenario(pdia)
{
    m_velocity = 0;
    m_dt = 1;
}

//-----------------------------------------------------------------------------
// Build the uniaxial loading scenario
// Cube with uniaxial velocity prescribed along x on left face and dilatation
// fixed on right face. Steady-state analysis.
bool FEFluidTangentUniaxialSS::Init()
{
    int i;
    vec3d r[8] = {
        vec3d(0,0,0), vec3d(1,0,0), vec3d(1,1,0), vec3d(0,1,0),
        vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1)
    };
    
/*    int BC[8][4] = {
        { 0,-1,-1, 0},{ 0,-1,-1,-1},{ 0,-1,-1,-1}, { 0,-1,-1, 0},
        { 0,-1,-1, 0},{ 0,-1,-1,-1},{ 0,-1,-1,-1}, { 0,-1,-1, 0}
    };*/
    int BC[8][4] = {
        { 0, 0, 0, 0},{ 0, 0, 0, 0},{ 0, 0, 0, 0}, { 0, 0, 0, 0},
        { 0, 0, 0, 0},{ 0, 0, 0, 0},{ 0, 0, 0, 0}, { 0, 0, 0, 0}
    };
    
    FEModel& fem = GetDiagnostic()->GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    pstep->m_nanalysis = FE_STEADY_STATE;
    
    int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
    const int dof_WX = fem.GetDOFIndex("wx");
    const int dof_WY = fem.GetDOFIndex("wy");
    const int dof_WZ = fem.GetDOFIndex("wz");
    const int dof_EF  = fem.GetDOFIndex("ef");
    
    // --- create the FE problem ---
    // create the mesh
    FEMesh& m = fem.GetMesh();
    m.CreateNodes(8);
    m.SetDOFS(MAX_DOFS);
    for (i=0; i<8; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_r0 = r[i];
        n.m_rid = -1;
        
        // set displacement BC's
        if (BC[i][0] == -1) fem.AddFixedBC(i, dof_WX);
        if (BC[i][1] == -1) fem.AddFixedBC(i, dof_WY);
        if (BC[i][2] == -1) fem.AddFixedBC(i, dof_WZ);
        if (BC[i][3] == -1) fem.AddFixedBC(i, dof_EF);
    }
    
    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    
    // create a fluid domain
    FEFluidDomain3D* pd = new FEFluidDomain3D(&fem);
    pd->SetMaterial(pmat);
    pd->Create(1, FE_HEX8G8);
    pd->SetMatID(0);
    m.AddDomain(pd);
    FESolidElement& el = pd->Element(0);
    el.SetID(1);
    for (i=0; i<8; ++i) el.m_node[i] = i;
    
    pd->CreateMaterialPointData();
    
    // Add a loadcurve
/*  FELoadCurve* plc = new FELoadCurve(new FELinearFunction(&fem, 1.0, 0.0));
    fem.AddLoadCurve(plc);
    
    // Add a prescribed BC
    int nd[4] = {0, 3, 4, 7};
    FEPrescribedDOF* pdc = new FEPrescribedDOF(&fem);
    fem.AddPrescribedBC(pdc);
    pdc->SetDOF(dof_WX).SetScale(m_velocity, 0);
    for (i = 0; i<4; ++i) pdc->AddNode(nd[i]);*/
    
    return true;
}

//-----------------------------------------------------------------------------
// Constructor
FEFluidTangentDiagnostic::FEFluidTangentDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
    m_pscn = 0;
    
    FEAnalysis* pstep = new FEAnalysis(&fem);
    
    // create a new solver
    FESolver* pnew_solver = fecore_new<FESolver>("fluid", &fem);
    assert(pnew_solver);
    pnew_solver->m_bsymm = false;
    pstep->SetFESolver(pnew_solver);
    
    fem.AddStep(pstep);
    fem.SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEDiagnosticScenario* FEFluidTangentDiagnostic::CreateScenario(const std::string& sname)
{
    if (sname == "fluid uni-axial") return (m_pscn = new FEFluidTangentUniaxial(this));
    else if (sname == "fluid uni-axial-SS") return (m_pscn = new FEFluidTangentUniaxialSS(this));
    return 0;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEFluidTangentDiagnostic::Init()
{
    if (m_pscn == 0) return false;
    
    if (m_pscn->Init() == false) return false;
    
    return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
// Helper function to print a matrix
void FEFluidTangentDiagnostic::print_matrix(matrix& m)
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
// Run the tangent diagnostic. After we run the FE model, we calculate
// the element stiffness matrix and compare that to a finite difference
// of the element residual.
bool FEFluidTangentDiagnostic::Run()
{
	Logfile::MODE oldmode = felog.SetMode(Logfile::LOG_FILE);
    
    // solve the problem
	felog.SetMode(Logfile::LOG_NEVER);
    FEModel& fem = GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
	fem.GetTime().timeIncrement = pstep->m_dt0 = dt;
    pstep->m_tstart = 0;
    pstep->m_tend = dt;
    pstep->m_final_time = dt;
    pstep->Activate();
    fem.Solve();
	felog.SetMode(Logfile::LOG_FILE);
    FETimeInfo tp;
    tp.timeIncrement = dt;
    tp.alpha = 1;
    tp.beta = 0.25;
    tp.gamma = 0.5;
    tp.alphaf = 1;
    tp.alpham = 1;
    
    FEMesh& mesh = fem.GetMesh();
    FEFluidDomain3D& bd = static_cast<FEFluidDomain3D&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = bd.Element(0);
    int N = mesh.Nodes();
    
    // set up the element stiffness matrix
    matrix k0(4*N, 4*N);
    k0.zero();
    bd.ElementStiffness(el, k0, tp);
    bd.ElementMassMatrix(el,k0, tp);
    
    // print the element stiffness matrix
    felog.printf("\nActual stiffness matrix:\n");
    print_matrix(k0);
    
    // now calculate the derivative of the residual
    matrix k1;
    deriv_residual(k1);
    
    // print the approximate element stiffness matrix
    felog.printf("\nApproximate stiffness matrix:\n");
    print_matrix(k1);
    
    // finally calculate the difference matrix
    felog.printf("\n");
    matrix kd(4*N, 4*N);
    double kmax = 0, kij;
    int i0 = -1, j0 = -1, i, j;
    for (i=0; i<4*N; ++i)
        for (j=0; j<4*N; ++j)
        {
            kd[i][j] = k0[i][j] - k1[i][j];
            kij = 100.0*fabs(kd[i][j] / k0[0][0]);
            if (kij > kmax)
            {
                kmax = kij;
                i0 = i;
                j0 = j;
            }
        }
    
    // print the difference
    felog.printf("\ndifference matrix:\n");
    print_matrix(kd);
    
    felog.SetMode(oldmode);
    
    felog.printf("\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);
    
    return (kmax < 1e-4);
}

//-----------------------------------------------------------------------------
// Calculate a finite difference approximation of the derivative of the
// element residual.
void FEFluidTangentDiagnostic::deriv_residual(matrix& ke)
{
    int i, j, nj;
    
    // get the solver
    FEModel& fem = GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
	fem.GetTime().timeIncrement = pstep->m_dt0 = dt;
    pstep->m_tstart = 0;
    pstep->m_tend = dt;
    pstep->m_final_time = dt;
    FEFluidSolver& solver = static_cast<FEFluidSolver&>(*pstep->GetFESolver());
    FETimeInfo tp;
    tp.timeIncrement = dt;
    tp.alpha = 1;
    tp.beta = 0.25;
    tp.gamma = 0.5;
    tp.alphaf = 1;
    tp.alpham = 1;

	// get the dof indices
	const int dof_WX = fem.GetDOFIndex("wx");
	const int dof_WY = fem.GetDOFIndex("wy");
	const int dof_WZ = fem.GetDOFIndex("wz");
	const int dof_EF  = fem.GetDOFIndex("ef");
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    FEFluidDomain3D& bd = static_cast<FEFluidDomain3D&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = bd.Element(0);
    int N = mesh.Nodes();
    
    // first calculate the initial residual
    vector<double> f0(4*N);
    zero(f0);
    bd.ElementInternalForce(el, f0, tp);
    bd.ElementInertialForce(el, f0, tp);
    
    // now calculate the perturbed residuals
    ke.resize(4*N, 4*N);
    ke.zero();
    double dx = 1e-8;
    vector<double> f1(4*N);
    for (j=0; j<4*N; ++j)
    {
        FENode& node = mesh.Node(el.m_node[j/4]);
        nj = j%4;
        
        switch (nj)
        {
            case 0: node.inc(dof_WX, dx); break;
            case 1: node.inc(dof_WY, dx); break;
            case 2: node.inc(dof_WZ, dx); break;
            case 3: node.inc(dof_EF, dx); break;
        }
        
        
		solver.UpdateModel();
        
        zero(f1);
        bd.ElementInternalForce(el, f1, tp);
        bd.ElementInertialForce(el, f1, tp);
        
        switch (nj)
        {
            case 0: node.dec(dof_WX, dx); break;
            case 1: node.dec(dof_WY, dx); break;
            case 2: node.dec(dof_WZ, dx); break;
            case 3: node.dec(dof_EF, dx); break;
        }
        
		solver.UpdateModel();
        
        for (i=0; i<4*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
