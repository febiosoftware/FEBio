// FEContactDiagnostic.cpp: implementation of the FEContactDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

#include "FEContactDiagnostic.h"
#include "FEBioMech/FENeoHookean.h"
#include "FEBioMech/FESolidSolver2.h"
#include "FEBioMech/FESlidingInterface.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FEBioMech/FEResidualVector.h"
#include "FECore/log.h"

void FEContactDiagnostic::print_matrix(DenseMatrix& m)
{
	int i, j;
	int N = m.Rows();
	int M = m.Columns();

	feLog("\n    ");
	for (i=0; i<N; ++i) feLog("%15d ", i);
	feLog("\n----");
	for (i=0; i<N; ++i) feLog("----------------", i);

	for (i=0; i<N; ++i)
	{
		feLog("\n%2d: ", i);
		for (j=0; j<M; ++j)
		{
			feLog("%15lg ", m(i,j));
		}
	}
	feLog("\n");
}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactDiagnostic::FEContactDiagnostic(FEModel& fem) : FEDiagnostic(fem)
{
	FEAnalysis* pstep = new FEAnalysis(&fem);
    fem.AddStep(pstep);
    fem.SetCurrentStep(pstep);
}

FEContactDiagnostic::~FEContactDiagnostic()
{

}

bool FEContactDiagnostic::Run()
{
	// get the solver
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();
	FESolidSolver2& solver = static_cast<FESolidSolver2&>(*pstep->GetFESolver());
	solver.Init();

	// make sure contact data is up to data
	fem.Update();

	// create the stiffness matrix
	solver.CreateStiffness(true);

	// get the stiffness matrix
	FEGlobalMatrix& K = *solver.GetStiffnessMatrix();
	SparseMatrix *pA = (SparseMatrix*)(&K);
	DenseMatrix& K0 = static_cast<DenseMatrix&>(*pA);

	// build the stiffness matrix
	K0.Zero();
	solver.ContactStiffness();
//	solver.StiffnessMatrix();

	print_matrix(K0);

	// calculate the derivative of the residual
	DenseMatrix K1;
	deriv_residual(K1);

	print_matrix(K1);

	// calculate difference matrix
	const int N = 48;
	DenseMatrix Kd; Kd.Create(N, N);
	double kij, kmax = 0, k0;
	int i0=-1, j0=-1;
	for (int i=0; i<N; ++i)
		for (int j=0; j<N; ++j)
		{
			Kd(i,j) = K1(i,j) - K0(i,j);
			k0 = fabs(K0(i,j));
			if (k0 < 1e-15) k0 = 1;
			kij = 100.0*fabs(Kd(i,j))/k0;
			if (kij > kmax) 
			{
				kmax = kij;
				i0 = i;
				j0 = j;
			}
		}

	print_matrix(Kd);

	feLog("\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);

	return (kmax < 1e-3);
}

//-----------------------------------------------------------------------------
bool FEContactDiagnostic::Init()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// --- create the geometry ---
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();

	// currently we simply assume a two-element contact problem
	// so we create two elements
	const double eps = 0.5;
	mesh.CreateNodes(16);
	mesh.SetDOFS(MAX_DOFS);
	mesh.Node( 0).m_r0 = vec3d(0,0,0);
	mesh.Node( 1).m_r0 = vec3d(1,0,0);
	mesh.Node( 2).m_r0 = vec3d(1,1,0);
	mesh.Node( 3).m_r0 = vec3d(0,1,0);
	mesh.Node( 4).m_r0 = vec3d(0,0,1);
	mesh.Node( 5).m_r0 = vec3d(1,0,1);
	mesh.Node( 6).m_r0 = vec3d(1,1,1);
	mesh.Node( 7).m_r0 = vec3d(0,1,1);
	mesh.Node( 8).m_r0 = vec3d(0,0,1-eps);
	mesh.Node( 9).m_r0 = vec3d(1,0,1-eps);
	mesh.Node(10).m_r0 = vec3d(1,1,1-eps);
	mesh.Node(11).m_r0 = vec3d(0,1,1-eps);
	mesh.Node(12).m_r0 = vec3d(0,0,2-eps);
	mesh.Node(13).m_r0 = vec3d(1,0,2-eps);
	mesh.Node(14).m_r0 = vec3d(1,1,2-eps);
	mesh.Node(15).m_r0 = vec3d(0,1,2-eps);

	for (int i=0; i<16; ++i)
	{
		FENode& node = mesh.Node(i);
		node.m_rt = node.m_r0;
	}

	// --- create a material ---
	FENeoHookean* pm = new FENeoHookean(&fem);
	pm->m_E = 1.0;
	pm->m_v = 0.45;
	fem.AddMaterial(pm);

	// get the one-and-only domain
	FEElasticSolidDomain* pbd = new FEElasticSolidDomain(&fem);
	pbd->SetMaterial(pm);
	pbd->Create(2, FE_HEX8G8);
	pbd->SetMatID(0);
	mesh.AddDomain(pbd);

	FESolidElement& el0 = pbd->Element(0);
	FESolidElement& el1 = pbd->Element(1);

	el0.SetID(1);
	el0.m_node[0] = 0;
	el0.m_node[1] = 1;
	el0.m_node[2] = 2;
	el0.m_node[3] = 3;
	el0.m_node[4] = 4;
	el0.m_node[5] = 5;
	el0.m_node[6] = 6;
	el0.m_node[7] = 7;

	el1.SetID(2);
	el1.m_node[0] = 8;
	el1.m_node[1] = 9;
	el1.m_node[2] = 10;
	el1.m_node[3] = 11;
	el1.m_node[4] = 12;
	el1.m_node[5] = 13;
	el1.m_node[6] = 14;
	el1.m_node[7] = 15;

	// --- create the sliding interface ---
	FESlidingInterface* ps = new FESlidingInterface(&fem);
	ps->m_atol = 0.1;
	ps->m_eps = 1;
	ps->m_btwo_pass = false;
	ps->m_nsegup = 0;
	FESlidingSurface& ms = ps->m_ms;
	ms.Create(1, FE_QUAD4NI);
	ms.Element(0).m_node[0] = 4;
	ms.Element(0).m_node[1] = 5;
	ms.Element(0).m_node[2] = 6;
	ms.Element(0).m_node[3] = 7;
	FESlidingSurface& ss = ps->m_ss;
	ss.Create(1, FE_QUAD4NI);
	ss.Element(0).m_node[0] = 11;
	ss.Element(0).m_node[1] = 10;
	ss.Element(0).m_node[2] = 9;
	ss.Element(0).m_node[3] = 8;
	fem.AddSurfacePairConstraint(ps);

	// --- set fem data ---
	// Make sure we are using the LU solver
	FECoreKernel::GetInstance().SetDefaultSolver("LU");

	return FEDiagnostic::Init();
}

//-----------------------------------------------------------------------------
void FEContactDiagnostic::deriv_residual(DenseMatrix& K)
{
	// get the solver
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();
	FESolidSolver2& solver = static_cast<FESolidSolver2&>(*pstep->GetFESolver());

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	fem.Update();

	// first calculate the initial residual
	vector<double> R0; R0.assign(48, 0);
	vector<double> dummy(R0);
	FEResidualVector RHS0(fem, R0, dummy);
	solver.ContactForces(RHS0);
//	solver.Residual(RHS);

	// now calculate the perturbed residuals
	K.Create(48, 48);
	int i, j, nj;
	int N = mesh.Nodes();
	double dx = 1e-8;
	vector<double> R1(48);
	for (j=0; j<3*N; ++j)
	{
		FENode& node = mesh.Node(j/3);
		nj = j%3;

		switch (nj)
		{
		case 0: node.m_rt.x += dx; break;
		case 1: node.m_rt.y += dx; break;
		case 2: node.m_rt.z += dx; break;
		}

		fem.Update();

		zero(R1);
		FEResidualVector RHS1(fem, R1, dummy);
		solver.ContactForces(RHS1);
//		solver.Residual(R1);

		switch (nj)
		{
		case 0: node.m_rt.x -= dx; break;
		case 1: node.m_rt.y -= dx; break;
		case 2: node.m_rt.z -= dx; break;
		}

		fem.Update();

		for (i=0; i<3*N; ++i) K(i,j) = (1-(R1[i] - R0[i])/dx)-1;
	}
}
