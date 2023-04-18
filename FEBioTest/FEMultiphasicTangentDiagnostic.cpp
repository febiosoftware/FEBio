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
#include "FEMultiphasicTangentDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEBioMix/FEMultiphasicSolver.h"
#include "FEBioMix/FEMultiphasicSolidDomain.h"
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEFixedBC.h>
#include "FECore/FEInitialCondition.h"
#include <FECore/FELoadCurve.h>
#include "FECore/log.h"


//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEMultiphasicTangentUniaxial, FEDiagnosticScenario)
	ADD_PARAMETER(m_strain       , "solid_strain"  );
	ADD_PARAMETER(m_pressure     , "fluid_pressure");
	ADD_PARAMETER(m_dt           , "time_step"     );
	ADD_PARAMETER(m_concentration, "solute_concentration");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMultiphasicTangentUniaxial::FEMultiphasicTangentUniaxial(FEDiagnostic* pdia) : FEMultiphasicScenario(pdia)
{ 
	m_strain = 0;
    m_pressure = 0;
    m_concentration = 0;
}

//-----------------------------------------------------------------------------
// Build the uniaxial loading scenario
// Cube with face displaced along x, with zero fluid pressure on displaced face
// and impermeable fixed face.
bool FEMultiphasicTangentUniaxial::Init()
{
    int i, isol;
    vec3d r[8] = {
        vec3d(0,0,0), vec3d(1,0,0), vec3d(1,1,0), vec3d(0,1,0),
        vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1)
    };
    
    int BC[8][3] = {
        {-1,-1,-1},{ 0,-1,-1},{ 0, 0,-1}, {-1, 0,-1},
        {-1,-1, 0},{ 0,-1, 0},{ 0, 0, 0}, {-1, 0, 0}
    };
    
    // get the material
	FEModel& fem = *GetDiagnostic()->GetFEModel();

    FEMaterial* pmat = fem.GetMaterial(0);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
    assert(pmp);
    int nsol = pmp->Solutes();
    double osm = nsol*m_concentration;
    double pe = -pmp->m_Rgas*pmp->m_Tabs*osm;
    
    // --- create the FE problem ---
	// get the degrees of freedom
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	const int dof_x = fem.GetDOFIndex("x");
	const int dof_y = fem.GetDOFIndex("y");
	const int dof_z = fem.GetDOFIndex("z");
	const int dof_p = fem.GetDOFIndex("p");
	const int dof_c = fem.GetDOFIndex("concentration", 0);

    // create the mesh
    FEMesh& m = fem.GetMesh();
    m.CreateNodes(8);
	m.SetDOFS(MAX_DOFS);

	// create a node set
	FENodeSet* iset = new FENodeSet(&fem);
	for (int i = 0; i < 8; ++i) iset->Add(i);

	// add initial conditions
	for (isol = 0; isol<nsol; ++isol) {

		FEInitialDOF* pic = new FEInitialDOF(&fem, dof_c + isol, iset);
		pic->SetValue(m_concentration);

		fem.AddInitialCondition(pic);
	}

	FEInitialDOF* pip = new FEInitialDOF(&fem, dof_p, iset);
	pip->SetValue(pe);
	fem.AddInitialCondition(pip);

	// add boundary conditions
	FENodeSet* nset[3];
	nset[0] = new FENodeSet(&fem);
	nset[1] = new FENodeSet(&fem);
	nset[2] = new FENodeSet(&fem);
	for (i=0; i<8; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_rp = n.m_r0 = r[i];
        
        // set displacement BC's
        if (BC[i][0] == -1) nset[0]->Add(i);
        if (BC[i][1] == -1) nset[1]->Add(i);
        if (BC[i][2] == -1) nset[2]->Add(i);
    }

	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_x, nset[0]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_y, nset[1]));
	fem.AddBoundaryCondition(new FEFixedBC(&fem, dof_z, nset[2]));

    // create a multiphasic domain
    FEMultiphasicSolidDomain* pd = new FEMultiphasicSolidDomain(&fem);
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

    // Add a prescribed displacement BC along X
	FENodeSet* dc = new FENodeSet(&fem);
	dc->Add({ 1, 2, 5, 6 });

	FEPrescribedDOF* pdc = new FEPrescribedDOF(&fem, dof_x, dc);
	pdc->SetScale(d, 0);
	fem.AddBoundaryCondition(pdc);
	
    // Add a prescribed fluid pressure BC
	FEPrescribedDOF* ppc = new FEPrescribedDOF(&fem, dof_p, dc);
	ppc->SetScale(pe, 0);
	fem.AddBoundaryCondition(ppc);

    // Add prescribed solute concentration BC
    for (i=0; i<nsol; ++i) {
		FEPrescribedDOF* psc = new FEPrescribedDOF(&fem, dof_c + i, dc);
		psc->SetScale(m_concentration, 0);
		fem.AddBoundaryCondition(psc);
	}

	return true;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Constructor
FEMultiphasicTangentDiagnostic::FEMultiphasicTangentDiagnostic(FEModel* fem) : FEDiagnostic(fem)
{
	// make sure the correct module is active
	fem->SetActiveModule("multiphasic");

    m_pscn = 0;

	FEAnalysis* pstep = new FEAnalysis(fem);

	// create a new solver
	FESolver* pnew_solver = fecore_new<FESolver>("multiphasic", fem);
	assert(pnew_solver);
	pnew_solver->m_msymm = REAL_UNSYMMETRIC;
	pstep->SetFESolver(pnew_solver);

	fem->AddStep(pstep);
	fem->SetCurrentStep(pstep);
}

//-----------------------------------------------------------------------------
FEDiagnosticScenario* FEMultiphasicTangentDiagnostic::CreateScenario(const std::string& sname)
{
	if (sname == "multiphasic uni-axial") return (m_pscn = new FEMultiphasicTangentUniaxial(this));
	return 0;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEMultiphasicTangentDiagnostic::Init()
{
	if (m_pscn == 0) return false;

	if (m_pscn->Init() == false) return false;

	return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
// Helper function to print a matrix
void FEMultiphasicTangentDiagnostic::print_matrix(matrix& m)
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
bool FEMultiphasicTangentDiagnostic::Run()
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
    if (fem.Solve() == false)
	{
		fem.UnBlockLog();
		return false;
	}
	fem.UnBlockLog();
    
    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
    assert(pmp);
    int nsol = pmp->Solutes();
    int ndpn = 4+nsol;
    
    FEMesh& mesh = fem.GetMesh();
    FEMultiphasicSolidDomain& md = static_cast<FEMultiphasicSolidDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = md.Element(0);
    int N = mesh.Nodes();
    
    // set up the element stiffness matrix
    matrix k0(ndpn*N, ndpn*N);
    k0.zero();
    md.ElementMultiphasicStiffness(el, k0, false);
    double k0max = 0;
    for (int i=0; i<ndpn*N; ++i)
        for (int j=0; j<ndpn*N; ++j)
            if (fabs(k0[i][j]) > k0max) k0max = fabs(k0[i][j]);
    
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
            kij = 100.0*fabs(kd[i][j] / k0max);
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
void FEMultiphasicTangentDiagnostic::deriv_residual(matrix& ke)
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
	FEMultiphasicSolver& solver = static_cast<FEMultiphasicSolver&>(*pstep->GetFESolver());
    
    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
    assert(pmp);
    int nsol = pmp->Solutes();
    int ndpn = 4+nsol;
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
    const int dof_x = fem.GetDOFIndex("x");
    const int dof_y = fem.GetDOFIndex("y");
    const int dof_z = fem.GetDOFIndex("z");
	const int dof_p = fem.GetDOFIndex("p");
	const int dof_c = fem.GetDOFIndex("concentration", 0);
    
    FEMultiphasicSolidDomain& md = static_cast<FEMultiphasicSolidDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = md.Element(0);
    int N = mesh.Nodes();
    
    // first calculate the initial residual
    vector<double> f0(ndpn*N);
    vector< vector<double> > f0c(nsol,vector<double>(N));
    zero(f0);
    md.ElementInternalForce(el, f0);
    
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
            case 0: node.add(dof_x, dx); node.m_rt.x += dx; break;
            case 1: node.add(dof_y, dx); node.m_rt.y += dx; break;
            case 2: node.add(dof_z, dx); node.m_rt.z += dx; break;
            case 3: node.add(dof_p, dx); break;
            default: node.add(dof_c + nj-4, dx); break;
        }
        
        
		fem.Update();
        
        zero(f1);
        md.ElementInternalForce(el, f1);
        
        switch (nj)
        {
            case 0: node.sub(dof_x, dx); node.m_rt.x -= dx; break;
            case 1: node.sub(dof_y, dx); node.m_rt.y -= dx; break;
            case 2: node.sub(dof_z, dx); node.m_rt.z -= dx; break;
            case 3: node.sub(dof_p, dx); break;
            default: node.sub(dof_c + nj-4, dx); break;
        }
        
		fem.Update();
        
        for (i=0; i<ndpn*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
