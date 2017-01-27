//
//  FEMultiphasicTangentDiagnostic.cpp
//  FEBio2
//
//  Created by Gerard Ateshian on 8/21/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEMultiphasicTangentDiagnostic.h"
#include "FETangentDiagnostic.h"
#include "FEBioLib/FEBox.h"
#include "FEBioMix/FEMultiphasicSolver.h"
#include "FEBioMix/FEMultiphasicDomain.h"
#include "FECore/BC.h"
#include "FECore/FEInitialCondition.h"
#include "FECore/LoadCurve.h"
#include "FECore/log.h"


//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEMultiphasicTangentUniaxial, FEDiagnosticScenario)
	ADD_PARAMETER(m_strain       , FE_PARAM_DOUBLE, "solid_strain"  );
	ADD_PARAMETER(m_pressure     , FE_PARAM_DOUBLE, "fluid_pressure");
	ADD_PARAMETER(m_dt           , FE_PARAM_DOUBLE, "time_step"     );
	ADD_PARAMETER(m_concentration, FE_PARAM_DOUBLE, "solute_concentration");
END_PARAMETER_LIST();

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
	FEModel& fem = GetDiagnostic()->GetFEModel();

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

	// add initial conditions
	for (isol = 0; isol<nsol; ++isol) {
		FEInitialBC* pic = new FEInitialBC(&fem);
		pic->SetDOF(dof_c + isol);
		for (i=0; i<8; ++i) pic->Add(i, m_concentration);
		fem.AddInitialCondition(pic);
	}

	FEInitialBC* pip = new FEInitialBC(&fem);
	pip->SetDOF(dof_p);
	for (i=0; i<8; ++i) pip->Add(i, pe);
	fem.AddInitialCondition(pip);

	// add boundary conditions
    for (i=0; i<8; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_rp = n.m_r0 = r[i];
        n.m_rid = -1;
        
        // set displacement BC's
        if (BC[i][0] == -1) fem.AddFixedBC(i, dof_x);
        if (BC[i][1] == -1) fem.AddFixedBC(i, dof_y);
        if (BC[i][2] == -1) fem.AddFixedBC(i, dof_z);
    }
    
    // create a multiphasic domain
    FEMultiphasicDomain* pd = new FEMultiphasicDomain(&fem);
	pd->SetMaterial(pmat);
    pd->Create(1, FE_HEX8G8);
	pd->SetMatID(0);
    m.AddDomain(pd);
    FESolidElement& el = pd->Element(0);
    el.SetID(1);
    for (i=0; i<8; ++i) el.m_node[i] = i;
    
    pd->CreateMaterialPointData();
    
    // convert strain to a displacement
    double d = sqrt(2*m_strain+1) - 1;
    
    // Add a loadcurve
    FELoadCurve* plc = new FELoadCurve;
    plc->Add(0, 0);
    plc->Add(1, 1);
    fem.AddLoadCurve(plc);
    
    // Add a prescribed displacement BC along X
    int nd[4] = {1, 2, 5, 6};
	FEPrescribedDOF* pdc = new FEPrescribedDOF(&fem);
    fem.AddPrescribedBC(pdc);
    pdc->SetDOF(dof_x).SetScale(d, 0);
    for (i = 0; i<4; ++i) pdc->AddNode(nd[i]);
    
    // Add a prescribed fluid pressure BC
	FEPrescribedDOF* ppc = new FEPrescribedDOF(&fem);
    fem.AddPrescribedBC(ppc);
    ppc->SetDOF(dof_p).SetScale(pe, 0);
    for (i = 0; i<4; ++i) ppc->AddNode(nd[i]);
    
    // Add prescribed solute concentration BC
    for (i=0; i<nsol; ++i) {
		FEPrescribedDOF* psc = new FEPrescribedDOF(&fem);
        fem.AddPrescribedBC(psc);
        psc->SetDOF(dof_c+i).SetScale(m_concentration, 0);
        for (i = 0; i<4; ++i) psc->AddNode(nd[i]);
    }

	return true;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Constructor
FEMultiphasicTangentDiagnostic::FEMultiphasicTangentDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
    m_pscn = 0;

	FEAnalysis* pstep = new FEAnalysis(&fem);

	// create a new solver
	FESolver* pnew_solver = fecore_new<FESolver>(FESOLVER_ID, "multiphasic", &fem);
	assert(pnew_solver);
	pnew_solver->m_bsymm = false;
	pstep->SetFESolver(pnew_solver);

	fem.AddStep(pstep);
	fem.m_nStep = 0;
	fem.SetCurrentStep(pstep);
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
bool FEMultiphasicTangentDiagnostic::Run()
{
	Logfile::MODE oldmode = felog.SetMode(Logfile::LOG_FILE);
    
    FEModel& fem = GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
    pstep->m_dt = pstep->m_dt0 =dt;
    pstep->m_tstart = 0;
    pstep->m_tend = dt;
    pstep->m_final_time = dt;

    // solve the problem
	felog.SetMode(Logfile::LOG_NEVER);
    if (fem.Solve() == false)
	{
		felog.SetMode(oldmode);
		return false;
	}
	felog.SetMode(Logfile::LOG_FILE);
    
    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
    assert(pmp);
    int nsol = pmp->Solutes();
    int ndpn = 4+nsol;
    
    FEMesh& mesh = fem.GetMesh();
    FEMultiphasicDomain& md = static_cast<FEMultiphasicDomain&>(mesh.Domain(0));
    
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
            kij = 100.0*fabs(kd[i][j] / k0max);
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
void FEMultiphasicTangentDiagnostic::deriv_residual(matrix& ke)
{
    int i, j, k, nj, isol;
    
    // get the solver
	FEModel& fem = GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
    double dt = m_pscn->m_dt;
    pstep->m_dt = pstep->m_dt0 =dt;
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
    
    FEMultiphasicDomain& md = static_cast<FEMultiphasicDomain&>(mesh.Domain(0));
    
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
            case 0: node.inc(dof_x, dx); node.m_rt.x += dx; break;
            case 1: node.inc(dof_y, dx); node.m_rt.y += dx; break;
            case 2: node.inc(dof_z, dx); node.m_rt.z += dx; break;
            case 3: node.inc(dof_p, dx); break;
            default: node.inc(dof_c + nj-4, dx); break;
        }
        
        
        solver.UpdateStresses();
        
        zero(f1);
        md.ElementInternalForce(el, f1);
        
        switch (nj)
        {
            case 0: node.dec(dof_x, dx); node.m_rt.x -= dx; break;
            case 1: node.dec(dof_y, dx); node.m_rt.y -= dx; break;
            case 2: node.dec(dof_z, dx); node.m_rt.z -= dx; break;
            case 3: node.dec(dof_p, dx); break;
            default: node.dec(dof_c + nj-4, dx); break;
        }
        
        solver.UpdateStresses();
        
        for (i=0; i<ndpn*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
