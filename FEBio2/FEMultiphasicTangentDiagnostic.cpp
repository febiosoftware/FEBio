//
//  FEMultiphasicTangentDiagnostic.cpp
//  FEBio2
//
//  Created by Gerard Ateshian on 8/21/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEMultiphasicTangentDiagnostic.h"
#include "stdafx.h"
#include "FETangentDiagnostic.h"
#include "FEBioLib/FEBox.h"
#include "FEBioMix/FEMultiphasicSolver.h"
#include "FEBioMix/FEMultiphasicDomain.h"
#include "FECore/log.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Constructor
FEMultiphasicTangentDiagnostic::FEMultiphasicTangentDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
    m_strain = 0;
    m_pressure = 0;
    m_concentration = 0;
    m_dt = 1;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEMultiphasicTangentDiagnostic::Init()
{
    switch (m_scn)
    {
        case TDS_MULTIPHASIC_UNIAXIAL: BuildUniaxial(); break;
        default:
            return false;
    }
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
// Build the uniaxial loading scenario
// Cube with face displaced along x, with prescribed fluid pressure
// and solute concentrations on displaced face
// and impermeable fixed face.
void FEMultiphasicTangentDiagnostic::BuildUniaxial()
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
    FEMaterial* pmat = m_fem.GetMaterial(0);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
    assert(pmp);
    int nsol = pmp->Solutes();
    double osm = nsol*m_concentration;
    double pe = -pmp->m_Rgas*pmp->m_Tabs*osm;
    
    // --- create the FE problem ---
    // create the mesh
    FEMesh& m = m_fem.GetMesh();
    m.CreateNodes(8);
    for (i=0; i<8; ++i)
    {
        FENode& n = m.Node(i);
        n.m_rt = n.m_rp = n.m_r0 = r[i];
        n.m_pt = n.m_p0 = pe;
        for (isol=0; isol<nsol; ++isol) {
            n.m_ct[isol] = n.m_cp[isol] = n.m_c0[isol] = m_concentration;
        }
        n.m_rid = -1;
        
        // set displacement BC's
        if (BC[i][0] == -1) m_fem.AddFixedBC(i, DOF_X);
        if (BC[i][1] == -1) m_fem.AddFixedBC(i, DOF_Y);
        if (BC[i][2] == -1) m_fem.AddFixedBC(i, DOF_Z);
    }
    
    // create a multiphasic domain
    FEMultiphasicDomain* pd = new FEMultiphasicDomain(&m, pmat);
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
    
    // Add a prescribed displacement BC along X
    int nd[4] = {1, 2, 5, 6};
    FEPrescribedBC* pdc = new FEPrescribedBC(&m_fem);
    m_fem.AddPrescribedBC(pdc);
    pdc->SetDOF(DOF_X).SetLoadCurveIndex(0).SetScale(d);
    for (i = 0; i<4; ++i) pdc->AddNode(nd[i]);
    
    // Add a prescribed fluid pressure BC
    FEPrescribedBC* ppc = new FEPrescribedBC(&m_fem);
    m_fem.AddPrescribedBC(ppc);
    ppc->SetDOF(DOF_P).SetLoadCurveIndex(0).SetScale(pe);
    for (i = 0; i<4; ++i) ppc->AddNode(nd[i]);
    
    // Add prescribed solute concentration BC
    for (i=0; i<nsol; ++i) {
        FEPrescribedBC* psc = new FEPrescribedBC(&m_fem);
        m_fem.AddPrescribedBC(psc);
        psc->SetDOF(DOF_C+i).SetLoadCurveIndex(0).SetScale(m_concentration);
        for (i = 0; i<4; ++i) psc->AddNode(nd[i]);
    }
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
    m_fem.Solve();
    felog.SetMode(Logfile::FILE_ONLY);
    
    // get the material
    FEMaterial* pmat = m_fem.GetMaterial(0);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
    assert(pmp);
    int nsol = pmp->Solutes();
    int ndpn = 4+nsol;
    
    FEMesh& mesh = m_fem.GetMesh();
    FEMultiphasicDomain& md = static_cast<FEMultiphasicDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = md.Element(0);
    int N = mesh.Nodes();
    
    // set up the element stiffness matrix
    matrix k0(ndpn*N, ndpn*N);
    k0.zero();
    md.ElementMultiphasicStiffness(el, k0, false, m_dt);
    
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
    FEAnalysis* pstep = m_fem.GetCurrentStep();
    FEMultiphasicSolver& solver = static_cast<FEMultiphasicSolver&>(*pstep->m_psolver);
    
    // get the material
    FEMaterial* pmat = m_fem.GetMaterial(0);
    FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pmat);
    assert(pmp);
    int nsol = pmp->Solutes();
    int ndpn = 4+nsol;
    
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    FEMultiphasicDomain& md = static_cast<FEMultiphasicDomain&>(mesh.Domain(0));
    
    // get the one and only element
    FESolidElement& el = md.Element(0);
    int N = mesh.Nodes();
    
    // first calculate the initial residual
    vector<double> f0(ndpn*N), f0u(3*N), f0p(N);
    vector< vector<double> > f0c(nsol,vector<double>(N));
    zero(f0u); zero(f0p); for (isol=0; isol<nsol; ++isol) zero(f0c[isol]);
    md.ElementInternalForce(el, f0u);
    md.ElementInternalFluidWork(el, f0p, m_dt);
    for (isol=0; isol<nsol; ++isol) md.ElementInternalSoluteWork(el, f0c[isol], m_dt, isol);
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
            case 3: node.m_pt   += dx; break;
            default: node.m_ct[nj-4] += dx; break;
        }
        
        
        solver.UpdateStresses();
        
        zero(f1u); zero(f1p); for (isol=0; isol<nsol; ++isol) zero(f1c[isol]);
        md.ElementInternalForce(el, f1u);
        md.ElementInternalFluidWork(el, f1p, m_dt);
        for (isol=0; isol<nsol; ++isol) md.ElementInternalSoluteWork(el, f1c[isol], m_dt, isol);
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
            case 3: node.m_pt   -= dx; break;
            default: node.m_ct[nj-4] -= dx; break;
        }
        
        solver.UpdateStresses();
        
        for (i=0; i<ndpn*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
    }
}
