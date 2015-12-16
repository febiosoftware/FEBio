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
    // create the mesh
    FEMesh& m = fem.GetMesh();
    m.CreateNodes(8);

	// get the degrees of freedom
	const int dof_x = fem.GetDOFIndex("x");
	const int dof_y = fem.GetDOFIndex("y");
	const int dof_z = fem.GetDOFIndex("z");
	const int dof_p = fem.GetDOFIndex("p");

	// add initial conditions
	for (isol = 0; isol<nsol; ++isol) {
		FEInitialBC* pic = new FEInitialBC(&fem);
		pic->SetDOF(DOF_C + isol);
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
    pd->create(1);
    m.AddDomain(pd);
    FESolidElement& el = pd->Element(0);
    el.SetType(FE_HEX8G8);
    el.m_nID = 1;
    el.SetMatID(0);
    for (i=0; i<8; ++i) el.m_node[i] = i;
    
    pd->InitMaterialPointData();
    
    // convert strain to a displacement
    double d = sqrt(2*m_strain+1) - 1;
    
    // Add a loadcurve
    FELoadCurve* plc = new FELoadCurve;
    plc->Add(0, 0);
    plc->Add(1, 1);
    fem.AddLoadCurve(plc);
    
    // Add a prescribed displacement BC along X
    int nd[4] = {1, 2, 5, 6};
    FEPrescribedBC* pdc = new FEPrescribedBC(&fem);
    fem.AddPrescribedBC(pdc);
    pdc->SetDOF(dof_x).SetLoadCurveIndex(0).SetScale(d);
    for (i = 0; i<4; ++i) pdc->AddNode(nd[i]);
    
    // Add a prescribed fluid pressure BC
    FEPrescribedBC* ppc = new FEPrescribedBC(&fem);
    fem.AddPrescribedBC(ppc);
    ppc->SetDOF(dof_p).SetLoadCurveIndex(0).SetScale(pe);
    for (i = 0; i<4; ++i) ppc->AddNode(nd[i]);
    
    // Add prescribed solute concentration BC
    for (i=0; i<nsol; ++i) {
        FEPrescribedBC* psc = new FEPrescribedBC(&fem);
        fem.AddPrescribedBC(psc);
        psc->SetDOF(DOF_C+i).SetLoadCurveIndex(0).SetScale(m_concentration);
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
    Logfile::MODE oldmode = felog.SetMode(Logfile::FILE_ONLY);
    
    // solve the problem
    felog.SetMode(Logfile::NEVER);
	FEModel& fem = GetFEModel();
    if (fem.Solve() == false)
	{
		felog.SetMode(oldmode);
		return false;
	}
    felog.SetMode(Logfile::FILE_ONLY);
    
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
	double dt = m_pscn->m_dt;
    md.ElementMultiphasicStiffness(el, k0, false, dt);
    
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
void FEMultiphasicTangentDiagnostic::deriv_residual(matrix& ke)
{
    int i, j, k, nj, isol;
    
    // get the solver
	FEModel& fem = GetFEModel();
    FEAnalysis* pstep = fem.GetCurrentStep();
	FEMultiphasicSolver& solver = static_cast<FEMultiphasicSolver&>(*pstep->GetFESolver());
    
    // get the material
    FEMaterial* pmat = fem.GetMaterial(0);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
    assert(pmp);
    int nsol = pmp->Solutes();
    int ndpn = 4+nsol;
    
    // get the mesh
    FEMesh& mesh = fem.GetMesh();
	const int dof_p = fem.GetDOFIndex("p");
    
    FEMultiphasicDomain& md = static_cast<FEMultiphasicDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = md.Element(0);
    int N = mesh.Nodes();
    
    // first calculate the initial residual
    vector<double> f0(ndpn*N), f0u(3*N), f0p(N);
    vector< vector<double> > f0c(nsol,vector<double>(N));
    zero(f0u); zero(f0p); for (isol=0; isol<nsol; ++isol) zero(f0c[isol]);
	double dt = m_pscn->m_dt;
    md.ElementInternalForce(el, f0u);
    md.ElementInternalFluidWork(el, f0p, dt);
    for (isol=0; isol<nsol; ++isol) md.ElementInternalSoluteWork(el, f0c[isol], dt, isol);
    for (i=0; i<N; ++i) {
        f0[ndpn*i  ] = f0u[3*i  ];
        f0[ndpn*i+1] = f0u[3*i+1];
        f0[ndpn*i+2] = f0u[3*i+2];
        f0[ndpn*i+3] = f0p[i    ];
        for (isol=0; isol<nsol; ++isol)
            f0[ndpn*i+4+isol] = f0c[isol][i];
    }
    
    // now calculate the perturbed residuals
    ke.resize(ndpn*N, ndpn*N);
    ke.zero();
    double dx = 1e-8;
    vector<double> f1(ndpn*N), f1u(3*N), f1p(N);
    vector< vector<double> > f1c(nsol,vector<double>(N));
    for (j=0; j<ndpn*N; ++j)
    {
        FENode& node = mesh.Node(el.m_node[j/ndpn]);
        nj = j%ndpn;
        
        switch (nj)
        {
            case 0: node.m_rt.x += dx; break;
            case 1: node.m_rt.y += dx; break;
            case 2: node.m_rt.z += dx; break;
            case 3: node.inc(dof_p, dx); break;
            default: node.inc(DOF_C + nj-4, dx); break;
        }
        
        
        solver.UpdateStresses();
        
        zero(f1u); zero(f1p); for (isol=0; isol<nsol; ++isol) zero(f1c[isol]);
        md.ElementInternalForce(el, f1u);
        md.ElementInternalFluidWork(el, f1p, dt);
        for (isol=0; isol<nsol; ++isol) md.ElementInternalSoluteWork(el, f1c[isol], dt, isol);
        for (k=0; k<N; ++k) {
            f1[ndpn*k  ] = f1u[3*k  ];
            f1[ndpn*k+1] = f1u[3*k+1];
            f1[ndpn*k+2] = f1u[3*k+2];
            f1[ndpn*k+3] = f1p[k    ];
            for (isol=0; isol<nsol; ++isol)
                f1[ndpn*k+4+isol] = f1c[isol][k];
        }
        
        switch (nj)
        {
            case 0: node.m_rt.x -= dx; break;
            case 1: node.m_rt.y -= dx; break;
            case 2: node.m_rt.z -= dx; break;
            case 3: node.dec(dof_p, dx); break;
            default: node.dec(DOF_C + nj-4, dx); break;
        }
        
        solver.UpdateStresses();
        
        for (i=0; i<ndpn*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
