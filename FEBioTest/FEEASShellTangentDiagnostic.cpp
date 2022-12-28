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
#include "FEEASShellTangentDiagnostic.h"
#include "FEBioMech/FESolidSolver2.h"
#include "FEBioMech/FEElasticEASShellDomain.h"
#include "FECore/log.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEEASShellTangentUnloaded, FEDiagnosticScenario)
	ADD_PARAMETER(m_strain, "strain");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
bool FEEASShellTangentUnloaded::Init()
{
    FEModel& fem = *GetDiagnostic()->GetFEModel();
    const int NELN = 4;
    
    int i;
    vec3d r[NELN] = {
        vec3d(0,0,0), vec3d(1,0,0), vec3d(1,1,0), vec3d(0,1,0)
    };
    
    vec3d D[NELN] = {
        vec3d(0,0,0.01), vec3d(0,0,0.01), vec3d(0,0,0.01), vec3d(0,0,0.01)
    };

    
    // get the degrees of freedom
    int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();

    // --- create the FE problem ---
    // create the mesh
    FEMesh& m = fem.GetMesh();
    m.CreateNodes(NELN);
    m.SetDOFS(MAX_DOFS);
    for (i=0; i<NELN; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_r0 = r[i];
        n.m_dt = n.m_d0 = D[i];
    }
    
    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    
    // create a solid domain
    FEElasticEASShellDomain* pd = new FEElasticEASShellDomain(&fem);
    pd->SetMaterial(pmat);
    pd->Create(1, FEElementLibrary::GetElementSpecFromType(FE_SHELL_QUAD4G8));
    pd->SetMatID(0);
    m.AddDomain(pd);
    FEShellElement& el = pd->Element(0);
    el.SetID(1);
    for (i=0; i<NELN; ++i) {
        el.m_node[i] = i;
        el.m_h0[i] = D[i].norm();
    }
    
    pd->CreateMaterialPointData();
    
    return true;
}

//-----------------------------------------------------------------------------
// Constructor
FEEASShellTangentDiagnostic::FEEASShellTangentDiagnostic(FEModel* fem) : FEDiagnostic(fem)
{
	// make sure the correct module is active
	fem->SetActiveModule("solid");

    m_pscn = 0;
    
    // create an analysis step
    FEAnalysis* pstep = new FEAnalysis(fem);
    
    // create a new solver
    FESolver* pnew_solver = fecore_new<FESolver>("solid", fem);
    assert(pnew_solver);
    pstep->SetFESolver(pnew_solver);
    
    // keep a pointer to the fem object
    fem->AddStep(pstep);
    fem->SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEDiagnosticScenario* FEEASShellTangentDiagnostic::CreateScenario(const std::string& sname)
{
    if (sname == "unloaded"   ) return (m_pscn = new FEEASShellTangentUnloaded   (this));
    return 0;
}

//-----------------------------------------------------------------------------
// Helper function to print a matrix
void FEEASShellTangentDiagnostic::print_matrix(matrix& m)
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
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEEASShellTangentDiagnostic::Init()
{
    // make sure we have a scenario
    if (m_pscn == 0) return false;
    
    // initialize the scenario
    if (m_pscn->Init() == false) return false;
    
    return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
// Run the tangent diagnostic. After we run the FE model, we calculate
// the element stiffness matrix and compare that to a finite difference
// of the element residual.
bool FEEASShellTangentDiagnostic::Run()
{
    // solve the problem
    FEModel& fem = *GetFEModel();
	fem.BlockLog();
    bool bret = fem.Solve();
	fem.UnBlockLog();
    if (bret == false) return false;
    
    FEMesh& mesh = fem.GetMesh();
    FEElasticEASShellDomain& bd = static_cast<FEElasticEASShellDomain&>(mesh.Domain(0));
    
    // set up the element stiffness matrix
    const int NELN = 4;
    const int NDPN = 6;
    const int NDOF = NELN*NDPN;
    matrix k0(NDOF, NDOF);
    k0.zero();
    bd.ElementStiffness(0, k0);
    
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
    matrix kd(NDOF, NDOF);
    double kmax = 0, kij;
    int i0 = -1, j0 = -1, i, j;
    for (i=0; i<NDOF; ++i)
        for (j=0; j<NDOF; ++j)
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
void FEEASShellTangentDiagnostic::deriv_residual(matrix& ke)
{
    // get the solver
    FEModel& fem = *GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    FESolidSolver2& solver = static_cast<FESolidSolver2&>(*pstep->GetFESolver());
    
    // get the degrees of freedom
    const int dof_X = fem.GetDOFIndex("x");
    const int dof_Y = fem.GetDOFIndex("y");
    const int dof_Z = fem.GetDOFIndex("z");
    const int dof_SX = fem.GetDOFIndex("sx");
    const int dof_SY = fem.GetDOFIndex("sy");
    const int dof_SZ = fem.GetDOFIndex("sz");

    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    
    FEElasticEASShellDomain& bd = static_cast<FEElasticEASShellDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FEShellElementNew& el = bd.ShellElement(0);
    
    // first calculate the initial residual
    const int NELN = 4;
    const int NDPN = 6;
    const int NDOF = NELN*NDPN;
    vector<double> f0(NDOF);
    zero(f0);
    bd.ElementInternalForce(el, f0);
    
    // now calculate the perturbed residuals
    ke.resize(NDOF, NDOF);
    ke.zero();
    int i, j, nj;
    int N = mesh.Nodes();
    double dx = 1e-8;
    vector<double> f1(NDOF);
    vector<double> ui(NDOF,0);
    for (j=0; j<NDPN*N; ++j)
    {
        FENode& node = mesh.Node(el.m_node[j/NDPN]);
        nj = j%NDPN;
        
        switch (nj)
        {
            case 0: node.add(dof_X, dx); node.m_rt.x += dx; break;
            case 1: node.add(dof_Y, dx); node.m_rt.y += dx; break;
            case 2: node.add(dof_Z, dx); node.m_rt.z += dx; break;
            case 3: node.add(dof_SX, dx); break;
            case 4: node.add(dof_SY, dx); break;
            case 5: node.add(dof_SZ, dx); break;
        }
        ui[j] += dx;
        
        solver.Update(ui);

        zero(f1);
        bd.ElementInternalForce(el, f1);
        
        switch (nj)
        {
            case 0: node.sub(dof_X, dx); node.m_rt.x -= dx; break;
            case 1: node.sub(dof_Y, dx); node.m_rt.y -= dx; break;
            case 2: node.sub(dof_Z, dx); node.m_rt.z -= dx; break;
            case 3: node.sub(dof_SX, dx); break;
            case 4: node.sub(dof_SY, dx); break;
            case 5: node.sub(dof_SZ, dx); break;
        }
        ui[j] -= dx;

        solver.Update(ui);

        for (i=0; i<NDPN*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
