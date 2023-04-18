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
#include "FEFluidFSITangentDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEBioLib/FEBox.h"
#include "FEBioFluid/FEFluidFSISolver.h"
#include "FEBioFluid/FEFluidFSIDomain3D.h"
#include "FEBioFluid/FEFluidFSIAnalysis.h"
#include "FECore/log.h"
#include <FECore/FEFixedBC.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidFSITangentUniaxial, FEFluidFSIScenario)
	ADD_PARAMETER(m_dilation, "fluid_dilation");
	ADD_PARAMETER(m_dt      , "time_step"     );
END_FECORE_CLASS();

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
    
    FEModel& fem = *GetDiagnostic()->GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    pstep->m_nanalysis = FEFluidFSIAnalysis::DYNAMIC;
//    pstep->m_nanalysis = FEFluidFSIAnalysis::STEADY_STATE;
    
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

	FENodeSet* nset[7] = { 0 };
	for (int i = 0; i < 7; ++i) nset[i] = new FENodeSet(&fem);
    for (i=0; i<8; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_r0 = r[i];
        
        // set displacement BC's
        if (BC[i][0] == -1) nset[0]->Add(i);
        if (BC[i][1] == -1) nset[1]->Add(i);
        if (BC[i][2] == -1) nset[2]->Add(i);
        if (BC[i][3] == -1) nset[3]->Add(i);
        if (BC[i][4] == -1) nset[4]->Add(i);
        if (BC[i][5] == -1) nset[5]->Add(i);
        if (BC[i][6] == -1) nset[6]->Add(i);
    }
    
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_X , nset[0]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_Y , nset[1]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_Z , nset[2]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_WX, nset[3]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_WY, nset[4]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_WZ, nset[5]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_EF, nset[6]));

    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    
    // create a fluid domain
    FEFluidFSIDomain3D* pd = new FEFluidFSIDomain3D(&fem);
    pd->SetMaterial(pmat);
    pd->Create(1, FEElementLibrary::GetElementSpecFromType(FE_HEX8G8));
    pd->SetMatID(0);
    m.AddDomain(pd);
    FESolidElement& el = pd->Element(0);
    el.SetID(1);
    for (i=0; i<8; ++i) el.m_node[i] = i;
    
    pd->CreateMaterialPointData();
    
    // Add a loadcurve
/*    FELoadCurve* plc = new FELoadCurve(new FELinearFunction(&fem, 1.0, 0.0));
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
FEFluidFSITangentDiagnostic::FEFluidFSITangentDiagnostic(FEModel* fem) : FEDiagnostic(fem)
{
    m_pscn = 0;

	// make sure the correct module is active
	fem->SetActiveModule("fluid-FSI");
    
    FEAnalysis* pstep = new FEAnalysis(fem);
    
    // create a new solver
    FESolver* pnew_solver = fecore_new<FESolver>("fluid-FSI", fem);
    assert(pnew_solver);
    pnew_solver->m_msymm = REAL_UNSYMMETRIC;
    FEFluidFSISolver* fluid_solver = dynamic_cast<FEFluidFSISolver*>(pnew_solver);
    fluid_solver->m_rhoi = 0;
    fluid_solver->m_pred = 0;
    pstep->SetFESolver(pnew_solver);
    
    fem->AddStep(pstep);
    fem->SetCurrentStep(pstep);
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
    
    feLog("\n    ");
    for (i=0; i<N; ++i) feLog("%15d ", i);
    feLog("\n----");
    for (i=0; i<N; ++i) feLog("----------------", i);
    
    for (i=0; i<N; ++i)
    {
        feLog("\n%2d: ", i);
        for (j=0; j<M; ++j)
        {
            feLog("%15lg ", m[i][j]);
        }
    }
    feLog("\n");
}

//-----------------------------------------------------------------------------
// Run the tangent diagnostic. After we run the FE model, we calculate
// the element stiffness matrix and compare that to a finite difference
// of the element residual.
bool FEFluidFSITangentDiagnostic::Run()
{
    // solve the problem
    FEModel& fem = *GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
	fem.GetTime().timeIncrement = pstep->m_dt0 = dt;
    pstep->m_tstart = 0;
    pstep->m_tend = dt;
    pstep->m_final_time = dt;
    pstep->Activate();
	fem.BlockLog();
    fem.Solve();
	fem.UnBlockLog();
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
    bd.ElementStiffness(el, k0);
    bd.ElementMassMatrix(el,k0);
    
    // print the element stiffness matrix
    feLog("\nActual stiffness matrix:\n");
    print_matrix(k0);
    
    // now calculate the derivative of the residual
    matrix k1;
    deriv_residual(k1);
    
    // print the approximate element stiffness matrix
    feLog("\nApproximate stiffness matrix:\n");
    print_matrix(k1);
    
    // finally calculate the difference matrix
    feLog("\n");
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
    feLog("\ndifference matrix:\n");
    print_matrix(kd);
    
    feLog("\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);
    
    return (kmax < 1e-4);
}

//-----------------------------------------------------------------------------
// Calculate a finite difference approximation of the derivative of the
// element residual.
void FEFluidFSITangentDiagnostic::deriv_residual(matrix& ke)
{
    int i, j, nj;
    
    // get the solver
    FEModel& fem = *GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
	fem.GetTime().timeIncrement = pstep->m_dt0 = dt;
    pstep->m_tstart = 0;
    pstep->m_tend = dt;
    pstep->m_final_time = dt;
    FEFluidFSISolver& solver = static_cast<FEFluidFSISolver&>(*pstep->GetFESolver());
    FETimeInfo& tp = fem.GetTime();
    tp.alpha = solver.m_alphaf;
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
    bd.ElementInternalForce(el, f0);
    bd.ElementInertialForce(el, f0);
    
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
            case 0: node.add(dof_X, dx); node.m_rt.x += dx; break;
            case 1: node.add(dof_Y, dx); node.m_rt.y += dx; break;
            case 2: node.add(dof_Z, dx); node.m_rt.z += dx; break;
            case 3: node.add(dof_WX, dx); break;
            case 4: node.add(dof_WY, dx); break;
            case 5: node.add(dof_WZ, dx); break;
            case 6: node.add(dof_EF, dx); break;
        }
        
        
		fem.Update();
        
        zero(f1);
        bd.ElementInternalForce(el, f1);
        bd.ElementInertialForce(el, f1);
        
        switch (nj)
        {
            case 0: node.sub(dof_X, dx); node.m_rt.x -= dx; break;
            case 1: node.sub(dof_Y, dx); node.m_rt.y -= dx; break;
            case 2: node.sub(dof_Z, dx); node.m_rt.z -= dx; break;
            case 3: node.sub(dof_WX, dx); break;
            case 4: node.sub(dof_WY, dx); break;
            case 5: node.sub(dof_WZ, dx); break;
            case 6: node.sub(dof_EF, dx); break;
        }
        
		fem.Update();
        
        for (i=0; i<ndpn*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
