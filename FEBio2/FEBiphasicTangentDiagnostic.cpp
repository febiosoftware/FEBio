//
//  FEBiphasicTangentDiagnostic.cpp
//  FEBio2
//
//  Created by Gerard Ateshian on 8/20/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEBiphasicTangentDiagnostic.h"
#include "stdafx.h"
#include "FETangentDiagnostic.h"
#include "FEBioLib/FEBox.h"
#include "FEBioMix/FEBiphasicSolver.h"
#include "FEBioMix/FEBiphasicSolidDomain.h"
#include "FECore/log.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Constructor
FEBiphasicTangentDiagnostic::FEBiphasicTangentDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
    m_strain = 0;
    m_pressure = 0;
    m_dt = 1;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEBiphasicTangentDiagnostic::Init()
{
    switch (m_scn)
    {
        case TDS_BIPHASIC_UNIAXIAL: BuildUniaxial(); break;
        default:
            return false;
    }
    return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
// Helper function to print a matrix
void FEBiphasicTangentDiagnostic::print_matrix(matrix& m)
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
// Build the uniaxial loading scenario
// Cube with face displaced along x, with zero fluid pressure on displaced face
// and impermeable fixed face.
void FEBiphasicTangentDiagnostic::BuildUniaxial()
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
    
    // --- create the FE problem ---
    // create the mesh
    FEMesh& m = m_fem.GetMesh();
    m.CreateNodes(8);
    for (i=0; i<8; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_r0 = r[i];
        n.m_pt = n.m_p0 = 0;
        n.m_rid = -1;
        
        // set displacement BC's
        if (BC[i][0] == -1) m_fem.AddFixedBC(i, DOF_X);
        if (BC[i][1] == -1) m_fem.AddFixedBC(i, DOF_Y);
        if (BC[i][2] == -1) m_fem.AddFixedBC(i, DOF_Z);
        if (BC[i][3] == -1) m_fem.AddFixedBC(i, DOF_P);
    }
    
    // get the material
    FEMaterial* pmat = m_fem.GetMaterial(0);
    
    // create a biphasic domain
    FEBiphasicSolidDomain* pd = new FEBiphasicSolidDomain(&m_fem);
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
    m_fem.AddLoadCurve(plc);
    
    // Add a prescribed BC
    int nd[4] = {1, 2, 5, 6};
    FEPrescribedBC* pdc = new FEPrescribedBC(&m_fem);
    m_fem.AddPrescribedBC(pdc);
    pdc->SetDOF(DOF_X).SetLoadCurveIndex(0).SetScale(d);
    for (i = 0; i<4; ++i) pdc->AddNode(nd[i]);
}

//-----------------------------------------------------------------------------
// Run the tangent diagnostic. After we run the FE model, we calculate
// the element stiffness matrix and compare that to a finite difference
// of the element residual.
bool FEBiphasicTangentDiagnostic::Run()
{
    Logfile::MODE oldmode = felog.SetMode(Logfile::FILE_ONLY);
    
    // solve the problem
    felog.SetMode(Logfile::NEVER);
    m_fem.Solve();
    felog.SetMode(Logfile::FILE_ONLY);
    
    FEMesh& mesh = m_fem.GetMesh();
    FEBiphasicSolidDomain& bd = static_cast<FEBiphasicSolidDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = bd.Element(0);
    int N = mesh.Nodes();
    
    // set up the element stiffness matrix
    matrix k0(4*N, 4*N);
    k0.zero();
    bd.ElementBiphasicStiffness(el, k0, false, m_dt);
    
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
void FEBiphasicTangentDiagnostic::deriv_residual(matrix& ke)
{
    int i, j, k, nj;
    
    // get the solver
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    FEBiphasicSolver& solver = static_cast<FEBiphasicSolver&>(*pstep->m_psolver);
    
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    FEBiphasicSolidDomain& bd = static_cast<FEBiphasicSolidDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = bd.Element(0);
    int N = mesh.Nodes();
    
    // first calculate the initial residual
    vector<double> f0(4*N), f0u(3*N), f0p(N);
    zero(f0u); zero(f0p);
    bd.ElementInternalForce(el, f0u);
    bd.ElementInternalFluidWork(el, f0p, m_dt);
    for (i=0; i<N; ++i) {
        f0[4*i  ] = f0u[3*i  ];
        f0[4*i+1] = f0u[3*i+1];
        f0[4*i+2] = f0u[3*i+2];
        f0[4*i+3] = f0p[i    ];
    }
    
    // now calculate the perturbed residuals
    ke.resize(4*N, 4*N);
    ke.zero();
    double dx = 1e-8;
    vector<double> f1(4*N), f1u(3*N), f1p(N);
    for (j=0; j<4*N; ++j)
    {
        FENode& node = mesh.Node(el.m_node[j/4]);
        nj = j%4;
        
        switch (nj)
        {
            case 0: node.m_rt.x += dx; break;
            case 1: node.m_rt.y += dx; break;
            case 2: node.m_rt.z += dx; break;
            case 3: node.m_pt   += dx; break;
        }
        
        
        solver.UpdateStresses();
        
        zero(f1u); zero(f1p);
        bd.ElementInternalForce(el, f1u);
        bd.ElementInternalFluidWork(el, f1p, m_dt);
        for (k=0; k<N; ++k) {
            f1[4*k  ] = f1u[3*k  ];
            f1[4*k+1] = f1u[3*k+1];
            f1[4*k+2] = f1u[3*k+2];
            f1[4*k+3] = f1p[k    ];
        }
        
        switch (nj)
        {
            case 0: node.m_rt.x -= dx; break;
            case 1: node.m_rt.y -= dx; break;
            case 2: node.m_rt.z -= dx; break;
            case 3: node.m_pt   -= dx; break;
        }
        
        solver.UpdateStresses();
        
        for (i=0; i<4*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
