//
//  FEFluidFSITangentDiagnostic.cpp
//  FEBioTest
//
//  Created by Gerard Ateshian on 8/16/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FEFluidFSITangentDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEBioLib/FEBox.h"
#include "FEBioFluid/FEFluidFSISolver.h"
#include "FEBioFluid/FEFluidFSIDomain3D.h"
#include "FECore/log.h"
#include <FECore/FEPrescribedBC.h>
#include "FECore/FEDataLoadCurve.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEFluidFSITangentUniaxial, FEFluidFSIScenario)
	ADD_PARAMETER(m_dilation, "fluid_dilation");
	ADD_PARAMETER(m_dt      , "time_step"     );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEFluidFSITangentUniaxial::FEFluidFSITangentUniaxial(FEDiagnostic* pdia) : FEFluidFSIScenario(pdia)
{
    m_dilation = 0;
    m_dt = 1;
}

//-----------------------------------------------------------------------------
// Build the uniaxial loading scenario
// Cube with uniaxial velocity prescribed along x on left face and dilatation
// fixed on right face. Dynamic analysis.
bool FEFluidFSITangentUniaxial::Init()
{
    int i;
    vec3d r[8] = {
        vec3d(-5.0000000e-01, -5.0000000e-01,  0.0000000e+00),
        vec3d( 5.0000000e-01, -5.0000000e-01,  0.0000000e+00),
        vec3d( 5.0000000e-01,  5.0000000e-01,  0.0000000e+00),
        vec3d(-5.0000000e-01,  5.0000000e-01,  0.0000000e+00),
        vec3d(-5.0000000e-01, -5.0000000e-01,  1.0000000e+00),
        vec3d( 5.0000000e-01, -5.0000000e-01,  1.0000000e+00),
        vec3d( 5.0000000e-01,  5.0000000e-01,  1.0000000e+00),
        vec3d(-5.0000000e-01,  5.0000000e-01,  1.0000000e+00)
    };
    
    int BC[8][7] = {
        { 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0}
    };
    
    FEModel& fem = GetDiagnostic()->GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    pstep->m_nanalysis = FE_DYNAMIC;
//    pstep->m_nanalysis = FE_STEADY_STATE;
    
    int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
    const int dof_X = fem.GetDOFIndex("x");
    const int dof_Y = fem.GetDOFIndex("y");
    const int dof_Z = fem.GetDOFIndex("z");
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
        if (BC[i][0] == -1) fem.AddFixedBC(i, dof_X);
        if (BC[i][1] == -1) fem.AddFixedBC(i, dof_Y);
        if (BC[i][2] == -1) fem.AddFixedBC(i, dof_Z);
        if (BC[i][3] == -1) fem.AddFixedBC(i, dof_WX);
        if (BC[i][4] == -1) fem.AddFixedBC(i, dof_WY);
        if (BC[i][5] == -1) fem.AddFixedBC(i, dof_WZ);
        if (BC[i][6] == -1) fem.AddFixedBC(i, dof_EF);
    }
    
    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    
    // create a fluid domain
    FEFluidFSIDomain3D* pd = new FEFluidFSIDomain3D(&fem);
    pd->SetMaterial(pmat);
    pd->Create(1, FE_HEX8G8);
    pd->SetMatID(0);
    m.AddDomain(pd);
    FESolidElement& el = pd->Element(0);
    el.SetID(1);
    for (i=0; i<8; ++i) el.m_node[i] = i;
    
    pd->CreateMaterialPointData();
    
    // Add a loadcurve
/*    FELoadCurve* plc = new FELinearRamp(1.0, 0.0);
    fem.AddLoadCurve(plc);
    
    // Add a prescribed BC
    int nd[4] = {0, 4, 3, 7};
    FEPrescribedDOF* pdc = new FEPrescribedDOF(&fem);
    fem.AddPrescribedBC(pdc);
    pdc->SetDOF(dof_EF).SetScale(m_dilation, 0);
    for (i = 0; i<4; ++i) pdc->AddNode(nd[i]);*/
    
    return true;
}

//-----------------------------------------------------------------------------
// Constructor
FEFluidFSITangentDiagnostic::FEFluidFSITangentDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
    m_pscn = 0;
    
    FEAnalysis* pstep = new FEAnalysis(&fem);
    
    // create a new solver
    FESolver* pnew_solver = fecore_new<FESolver>(FESOLVER_ID, "fluid-FSI", &fem);
    assert(pnew_solver);
    pnew_solver->m_bsymm = false;
    FEFluidFSISolver* fluid_solver = dynamic_cast<FEFluidFSISolver*>(pnew_solver);
    fluid_solver->m_rhoi = 0;
    fluid_solver->m_pred = 0;
    pstep->SetFESolver(pnew_solver);
    
    fem.AddStep(pstep);
    fem.SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEDiagnosticScenario* FEFluidFSITangentDiagnostic::CreateScenario(const std::string& sname)
{
    if (sname == "fluid uni-axial") return (m_pscn = new FEFluidFSITangentUniaxial(this));
    return 0;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEFluidFSITangentDiagnostic::Init()
{
    if (m_pscn == 0) return false;
    
    if (m_pscn->Init() == false) return false;
    
    return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
// Helper function to print a matrix
void FEFluidFSITangentDiagnostic::print_matrix(matrix& m)
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
bool FEFluidFSITangentDiagnostic::Run()
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
    FEFluidFSIDomain3D& bd = static_cast<FEFluidFSIDomain3D&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = bd.Element(0);
    int N = mesh.Nodes();
    
    // set up the element stiffness matrix
    const int ndpn = 7;
    matrix k0(ndpn*N, ndpn*N);
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
    matrix kd(ndpn*N, ndpn*N);
    double kmax = 0, kij;
    int i0 = -1, j0 = -1, i, j;
    for (i=0; i<ndpn*N; ++i)
        for (j=0; j<ndpn*N; ++j)
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
void FEFluidFSITangentDiagnostic::deriv_residual(matrix& ke)
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
    FEFluidFSISolver& solver = static_cast<FEFluidFSISolver&>(*pstep->GetFESolver());
    FETimeInfo& tp = fem.GetTime();
    tp.alpha = solver.m_alpha;
    tp.beta  = solver.m_beta;
    tp.gamma = solver.m_gamma;
    tp.alphaf = solver.m_alphaf;
    tp.alpham  = solver.m_alpham;
    
    // get the dof indices
    const int dof_X = fem.GetDOFIndex("x");
    const int dof_Y = fem.GetDOFIndex("y");
    const int dof_Z = fem.GetDOFIndex("z");
    const int dof_WX = fem.GetDOFIndex("wx");
    const int dof_WY = fem.GetDOFIndex("wy");
    const int dof_WZ = fem.GetDOFIndex("wz");
    const int dof_EF  = fem.GetDOFIndex("ef");
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    FEFluidFSIDomain3D& bd = static_cast<FEFluidFSIDomain3D&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = bd.Element(0);
    int N = mesh.Nodes();
    
    // first calculate the initial residual
    const int ndpn = 7;
    vector<double> f0(ndpn*N);
    zero(f0);
    bd.ElementInternalForce(el, f0, tp);
    bd.ElementInertialForce(el, f0, tp);
    
    // now calculate the perturbed residuals
    ke.resize(ndpn*N, ndpn*N);
    ke.zero();
    double dx = 1e-8;
    vector<double> f1(ndpn*N);
    for (j=0; j<ndpn*N; ++j)
    {
        FENode& node = mesh.Node(el.m_node[j/ndpn]);
        nj = j%ndpn;
        
        switch (nj)
        {
            case 0: node.inc(dof_X, dx); node.m_rt.x += dx; break;
            case 1: node.inc(dof_Y, dx); node.m_rt.y += dx; break;
            case 2: node.inc(dof_Z, dx); node.m_rt.z += dx; break;
            case 3: node.inc(dof_WX, dx); break;
            case 4: node.inc(dof_WY, dx); break;
            case 5: node.inc(dof_WZ, dx); break;
            case 6: node.inc(dof_EF, dx); break;
        }
        
        
		solver.UpdateModel();
        
        zero(f1);
        bd.ElementInternalForce(el, f1, tp);
        bd.ElementInertialForce(el, f1, tp);
        
        switch (nj)
        {
            case 0: node.dec(dof_X, dx); node.m_rt.x -= dx; break;
            case 1: node.dec(dof_Y, dx); node.m_rt.y -= dx; break;
            case 2: node.dec(dof_Z, dx); node.m_rt.z -= dx; break;
            case 3: node.dec(dof_WX, dx); break;
            case 4: node.dec(dof_WY, dx); break;
            case 5: node.dec(dof_WZ, dx); break;
            case 6: node.dec(dof_EF, dx); break;
        }
        
		solver.UpdateModel();
        
        for (i=0; i<ndpn*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
