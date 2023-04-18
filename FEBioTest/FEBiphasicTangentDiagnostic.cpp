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
#include "FEBiphasicTangentDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEBioMix/FEBiphasicSolver.h"
#include "FEBioMix/FEBiphasicSolidDomain.h"
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEFixedBC.h>
#include <FECore/FELoadCurve.h>
#include "FECore/log.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEBiphasicTangentUniaxial, FEBiphasicScenario)
	ADD_PARAMETER(m_strain  , "solid_strain"  );
	ADD_PARAMETER(m_pressure, "fluid_pressure");
	ADD_PARAMETER(m_dt      , "time_step"     );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEBiphasicTangentUniaxial::FEBiphasicTangentUniaxial(FEDiagnostic* pdia) : FEBiphasicScenario(pdia)
{ 
    m_strain = 0;
    m_pressure = 0;
}

//-----------------------------------------------------------------------------
// Build the uniaxial loading scenario
// Cube with face displaced along x, with zero fluid pressure on displaced face
// and impermeable fixed face.
bool FEBiphasicTangentUniaxial::Init()
{
    int i;
    vec3d r[8] = {
        vec3d(0,0,0), vec3d(1,0,0), vec3d(1,1,0), vec3d(0,1,0),
        vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1)
    };
    
    int BC[8][4] = {
        {-1,-1,-1, 0},{ 0,-1,-1,-1},{ 0, 0,-1,-1}, {-1, 0,-1, 0},
        {-1,-1, 0, 0},{ 0,-1, 0,-1},{ 0, 0, 0,-1}, {-1, 0, 0, 0}
    };
  
	FEModel& fem = *GetDiagnostic()->GetFEModel();
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	const int dof_x = fem.GetDOFIndex("x");
	const int dof_y = fem.GetDOFIndex("y");
	const int dof_z = fem.GetDOFIndex("z");
	const int dof_p = fem.GetDOFIndex("p");

    // --- create the FE problem ---
    // create the mesh
    FEMesh& m = fem.GetMesh();
    m.CreateNodes(8);
	m.SetDOFS(MAX_DOFS);

	FENodeSet* nset[4] = { 0 };
	for (int i = 0; i < 4; ++i) nset[i] = new FENodeSet(&fem);

    for (i=0; i<8; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_r0 = r[i];
        
        // set displacement BC's
        if (BC[i][0] == -1) nset[0]->Add(i);
        if (BC[i][1] == -1) nset[1]->Add(i);
        if (BC[i][2] == -1) nset[2]->Add(i);
        if (BC[i][3] == -1) nset[3]->Add(i);
    }
    
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_x, nset[0]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_y, nset[1]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_z, nset[2]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_p, nset[3]));

    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    
    // create a biphasic domain
    FEBiphasicSolidDomain* pd = new FEBiphasicSolidDomain(&fem);
	pd->SetMaterial(pmat);
    pd->Create(1, FEElementLibrary::GetElementSpecFromType(FE_HEX8G8));
	pd->SetMatID(0);
    m.AddDomain(pd);
    FESolidElement& el = pd->Element(0);
    el.SetID(1);
    for (i=0; i<8; ++i) el.m_node[i] = i;
    
    pd->CreateMaterialPointData();
    
    // convert strain to a displacement
    double d = sqrt(2*m_strain+1) - 1;
    
    // Add a loadcurve
	FELoadCurve* plc = new FELoadCurve(&fem);
	plc->Add(0.0, 0.0);
	plc->Add(1.0, 1.0);
    fem.AddLoadController(plc);
    
    // Add a prescribed BC
	FENodeSet* dc = new FENodeSet(&fem);
    dc->Add({1, 2, 5, 6});
    FEPrescribedDOF* pdc = new FEPrescribedDOF(&fem, dof_x, dc);
    pdc->SetScale(d, 0);
	fem.AddBoundaryCondition(pdc);

	return true;
}

//-----------------------------------------------------------------------------
// Constructor
FEBiphasicTangentDiagnostic::FEBiphasicTangentDiagnostic(FEModel* fem) : FEDiagnostic(fem)
{
	m_pscn = 0;

	// make sure the correct module is active
	fem->SetActiveModule("biphasic");

	FEAnalysis* pstep = new FEAnalysis(fem);

	// create a new solver
	FESolver* pnew_solver = fecore_new<FESolver>("biphasic", fem);
	assert(pnew_solver);
	pnew_solver->m_msymm = REAL_UNSYMMETRIC;
	pstep->SetFESolver(pnew_solver);

	fem->AddStep(pstep);
	fem->SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEDiagnosticScenario* FEBiphasicTangentDiagnostic::CreateScenario(const std::string& sname)
{
	if (sname == "biphasic uni-axial") return (m_pscn = new FEBiphasicTangentUniaxial(this));
	return 0;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEBiphasicTangentDiagnostic::Init()
{
	if (m_pscn == 0) return false;

	if (m_pscn->Init() == false) return false;
    
    return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
// Helper function to print a matrix
void FEBiphasicTangentDiagnostic::print_matrix(matrix& m)
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
bool FEBiphasicTangentDiagnostic::Run()
{
    FEModel& fem = *GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
	fem.GetTime().timeIncrement = pstep->m_dt0 = dt;
    pstep->m_tstart = 0;
    pstep->m_tend = dt;
    pstep->m_final_time = dt;
    
    // solve the problem
	fem.BlockLog();
    fem.Solve();
	fem.UnBlockLog();
    
    FEMesh& mesh = fem.GetMesh();
    FEBiphasicSolidDomain& bd = static_cast<FEBiphasicSolidDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = bd.Element(0);
    int N = mesh.Nodes();
    
    // set up the element stiffness matrix
    matrix k0(4*N, 4*N);
    k0.zero();
    bd.ElementBiphasicStiffness(el, k0, false);
    
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
    feLog("\ndifference matrix:\n");
    print_matrix(kd);
    
    feLog("\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);
    
    return (kmax < 1e-4);
}

//-----------------------------------------------------------------------------
// Calculate a finite difference approximation of the derivative of the
// element residual.
void FEBiphasicTangentDiagnostic::deriv_residual(matrix& ke)
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

	// get the DOFs
	const int dof_x = fem.GetDOFIndex("x");
	const int dof_y = fem.GetDOFIndex("y");
	const int dof_z = fem.GetDOFIndex("z");
	const int dof_p = fem.GetDOFIndex("p");
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    FEBiphasicSolidDomain& bd = static_cast<FEBiphasicSolidDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = bd.Element(0);
    int N = mesh.Nodes();
    
    // first calculate the initial residual
    vector<double> f0(4*N);
    zero(f0);
    bd.ElementInternalForce(el, f0);
    
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
            case 0: node.add(dof_x, dx); node.m_rt.x += dx; break;
            case 1: node.add(dof_y, dx); node.m_rt.y += dx; break;
            case 2: node.add(dof_z, dx); node.m_rt.z += dx; break;
            case 3: node.add(dof_p, dx); break;
        }
        
		fem.Update();
        
        zero(f1);
        bd.ElementInternalForce(el, f1);
        
        switch (nj)
        {
            case 0: node.sub(dof_x, dx); node.m_rt.x -= dx; break;
            case 1: node.sub(dof_y, dx); node.m_rt.y -= dx; break;
            case 2: node.sub(dof_z, dx); node.m_rt.z -= dx; break;
            case 3: node.sub(dof_p, dx); break;
        }
        
		fem.Update();
        
        for (i=0; i<4*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
