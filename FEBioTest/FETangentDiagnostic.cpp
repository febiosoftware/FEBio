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
#include "FETangentDiagnostic.h"
#include "FEBioMech/FESolidSolver2.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEFixedBC.h>
#include <FECore/FELoadCurve.h>
#include "FECore/log.h"
#include <FECore/FECoreKernel.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// Helper function to print a matrix
void print_matrix(FILE* fp, matrix& m)
{
	int i, j;
	int N = m.rows();
	int M = m.columns();

	fprintf(fp, "\n    ");
	for (i=0; i<N; ++i) fprintf(fp, "%15d ", i);
	fprintf(fp, "\n----");
	for (i=0; i<N; ++i) fprintf(fp, "----------------");

	for (i=0; i<N; ++i)
	{
		fprintf(fp, "\n%2d: ", i);
		for (j=0; j<M; ++j)
		{
			fprintf(fp, "%15lg ", m[i][j]);
		}
	}
	fprintf(fp, "\n");
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FETangentUniaxial, FEDiagnosticScenario)
	ADD_PARAMETER(m_strain, "strain");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
bool FETangentUniaxial::Init()
{
	FEModel& fem = *GetDiagnostic()->GetFEModel();

	int i;
	vec3d r[8] = {
		vec3d(0,0,0), vec3d(1,0,0), vec3d(1,1,0), vec3d(0,1,0),
		vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1)
	};

	int BC[8][3] = {
		{-1,-1,-1},{ 0,-1,-1},{ 0, 0,-1}, {-1, 0,-1},
		{-1,-1, 0},{ 0,-1, 0},{ 0, 0, 0}, {-1, 0, 0}
	};

	// get the degrees of freedom
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");

	// --- create the FE problem ---
	// create the mesh
	FEMesh& m = fem.GetMesh();
	m.CreateNodes(8);
	m.SetDOFS(MAX_DOFS);

	FENodeSet* nset[3];
	nset[0] = new FENodeSet(&fem);
	nset[1] = new FENodeSet(&fem);
	nset[2] = new FENodeSet(&fem);

	for (i=0; i<8; ++i)
	{
		FENode& n = m.Node(i);
		n.m_rt = n.m_r0 = r[i];

		// set displacement BC's
		if (BC[i][0] == -1) nset[0]->Add(i);
		if (BC[i][1] == -1) nset[1]->Add(i);
		if (BC[i][2] == -1) nset[2]->Add(i);
	}

	FEFixedBC* bc_x = new FEFixedBC(&fem, dof_X, nset[0]); fem.AddBoundaryCondition(bc_x);
	FEFixedBC* bc_y = new FEFixedBC(&fem, dof_Y, nset[1]); fem.AddBoundaryCondition(bc_y);
	FEFixedBC* bc_z = new FEFixedBC(&fem, dof_Z, nset[2]); fem.AddBoundaryCondition(bc_z);

	// get the material
	FEMaterial* pmat = fem.GetMaterial(0);

	FE_Element_Spec es;
	es.eclass = FE_Element_Class::FE_ELEM_SOLID;
	es.eshape = FE_Element_Shape::ET_HEX8;
	es.etype = FE_Element_Type::FE_HEX8G8;

	// create a solid domain
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	FEElasticSolidDomain* pd = dynamic_cast<FEElasticSolidDomain*>(fecore.CreateDomain(es, &m, pmat));
	pd->Create(1, es);
	pd->SetMatID(0);
	m.AddDomain(pd);
	FESolidElement& el = pd->Element(0);
	el.SetID(1);
	for (i=0; i<8; ++i) el.m_node[i] = i;

	pd->CreateMaterialPointData();

	// convert strain to a displacement
	double d = sqrt(2*m_strain+1) - 1;

	// get the simulation time
	FEAnalysis* step = fem.GetStep(0);
	double tend = step->m_ntime * step->m_dt0;

	// Add a loadcurve
	FELoadCurve* plc = new FELoadCurve(&fem);
	plc->Add(0.0, 0.0);
	plc->Add(tend, 1.0);
	fem.AddLoadController(plc);

	// Add a prescribed BC
	int nd[4] = {1, 2, 5, 6};
	FENodeSet* dc = new FENodeSet(&fem);
	for (i = 0; i<4; ++i) dc->Add(nd[i]);

	FEPrescribedDOF* pdc = new FEPrescribedDOF(&fem, dof_X, dc);
	pdc->SetScale(d, 0);
	fem.AddBoundaryCondition(pdc);

	return true;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FETangentSimpleShear, FEDiagnosticScenario)
	ADD_PARAMETER(m_strain, "strain");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
bool FETangentSimpleShear::Init()
{
	FEModel& fem = *GetDiagnostic()->GetFEModel();

	int i;
	vec3d r[8] = {
		vec3d(0,0,0), vec3d(1,0,0), vec3d(1,1,0), vec3d(0,1,0),
		vec3d(0,0,1), vec3d(1,0,1), vec3d(1,1,1), vec3d(0,1,1)
	};

	int BC[8][3] = {
		{-1,-1,-1},{-1,-1,-1},{-1, 0,-1}, {-1, 0,-1},
		{ 0,-1,-1},{ 0,-1,-1},{ 0, 0,-1}, { 0, 0,-1}
	};

	// get the degrees of freedom
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");

	// --- create the FE problem ---
	// create the mesh
	FEMesh& m = fem.GetMesh();
	m.CreateNodes(8);
	m.SetDOFS(MAX_DOFS);
	FENodeSet* nset[3];
	nset[0] = new FENodeSet(&fem);
	nset[1] = new FENodeSet(&fem);
	nset[2] = new FENodeSet(&fem);

	for (i = 0; i<8; ++i)
	{
		FENode& n = m.Node(i);
		n.m_rt = n.m_r0 = r[i];

		// set displacement BC's
		if (BC[i][0] == -1) nset[0]->Add(i);
		if (BC[i][1] == -1) nset[1]->Add(i);
		if (BC[i][2] == -1) nset[2]->Add(i);
	}

	FEFixedBC* bc_x = new FEFixedBC(&fem, dof_X, nset[0]); fem.AddBoundaryCondition(bc_x);
	FEFixedBC* bc_y = new FEFixedBC(&fem, dof_Y, nset[1]); fem.AddBoundaryCondition(bc_y);
	FEFixedBC* bc_z = new FEFixedBC(&fem, dof_Z, nset[2]); fem.AddBoundaryCondition(bc_z);

	// get the material
	FEMaterial* pmat = fem.GetMaterial(0);

	// create a solid domain
	FEElasticSolidDomain* pd = new FEElasticSolidDomain(&fem);
	pd->SetMaterial(pmat);
	pd->Create(1, FEElementLibrary::GetElementSpecFromType(FE_HEX8G8));
	pd->SetMatID(0);
	m.AddDomain(pd);
	FESolidElement& el = pd->Element(0);
	el.SetID(1);
	for (i=0; i<8; ++i) el.m_node[i] = i;

	pd->CreateMaterialPointData();

	// convert strain to a displacement
	double d = 2*m_strain;

	// Add a loadcurve
	FELoadCurve* plc = new FELoadCurve(&fem);
	plc->Add(0.0, 0.0);
	plc->Add(1.0, 1.0);
	fem.AddLoadController(plc);

	// Add a prescribed BC
	int nd[4] = { 4, 5, 6, 7 };
	FENodeSet* dc = new FENodeSet(&fem);
	for (i = 0; i<4; ++i) dc->Add(nd[i]);

	FEPrescribedDOF* pdc = new FEPrescribedDOF(&fem, dof_X, dc);
	pdc->SetScale(d, 0);
	fem.AddBoundaryCondition(pdc);

	return true;
}

//-----------------------------------------------------------------------------
// Constructor
FETangentDiagnostic::FETangentDiagnostic(FEModel* fem) : FEDiagnostic(fem)
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
FEDiagnosticScenario* FETangentDiagnostic::CreateScenario(const std::string& sname)
{
	if (sname == "uni-axial"   ) return (m_pscn = new FETangentUniaxial   (this));
	if (sname == "simple shear") return (m_pscn = new FETangentSimpleShear(this));
	return 0;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FETangentDiagnostic::Init()
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
bool FETangentDiagnostic::Run()
{
	// solve the problem
	FEModel& fem = *GetFEModel();
	fem.BlockLog();
	bool bret = fem.Solve();
	fem.UnBlockLog();
	if (bret == false)
	{
		feLogError("FEBio error terminated. Aborting diagnostic.\n");
		return false;
	}

	FEMesh& mesh = fem.GetMesh();
	FEElasticSolidDomain& bd = static_cast<FEElasticSolidDomain&>(mesh.Domain(0));

	// get the one and only element
	FESolidElement& el = bd.Element(0);

	// set up the element stiffness matrix
	matrix k0(24, 24);
	k0.zero();
	bd.ElementStiffness(fem.GetTime(), 0, k0);

	// create a file name for the tangent log file
	std::string febFile = GetFileName();
	string fileName = febFile;
	size_t n = fileName.rfind('.');
	if (n != string::npos) fileName.erase(n);
	fileName.append("_out.log");

	FILE* fp = fopen(fileName.c_str(), "wt");

	fprintf(fp, "FEBio Tangent Diagnostics Results:\n");
	fprintf(fp, "==================================\n");
	fprintf(fp, "Diagnostics file: %s\n\n", febFile.c_str());

	// print the element stiffness matrix
	fprintf(fp, "\nActual stiffness matrix:\n");
	print_matrix(fp, k0);

	// now calculate the derivative of the residual
	matrix k1;
	deriv_residual(k1);

	// print the approximate element stiffness matrix
	fprintf(fp, "\nApproximate stiffness matrix:\n");
	print_matrix(fp, k1);

	// finally calculate the difference matrix
	fprintf(fp, "\n");
	matrix kd(24, 24);
	double kmax = 0, kij;
	int i0 = -1, j0 = -1, i, j;
	for (i=0; i<24; ++i)
		for (j=0; j<24; ++j)
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
	fprintf(fp, "\ndifference matrix:\n");
	print_matrix(fp, kd);

	fprintf(fp, "\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);
	fprintf(stderr, "\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);

	fclose(fp);

	return (kmax < 1e-4);
}

//-----------------------------------------------------------------------------
// Calculate a finite difference approximation of the derivative of the
// element residual.
void FETangentDiagnostic::deriv_residual(matrix& ke)
{
	// get the solver
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();
	FESolidSolver2& solver = static_cast<FESolidSolver2&>(*pstep->GetFESolver());

	// get the degrees of freedom
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	FEElasticSolidDomain& bd = static_cast<FEElasticSolidDomain&>(mesh.Domain(0));

	// get the one and only element
	FESolidElement& el = bd.Element(0);

	// first calculate the initial residual
	vector<double> f0(24);
	zero(f0);
	bd.ElementInternalForce(el, f0);

	// now calculate the perturbed residuals
	ke.resize(24, 24);
	ke.zero();
	int i, j, nj;
	int N = mesh.Nodes();
	double dx = 1e-8;
	vector<double> f1(24);
	for (j=0; j<3*N; ++j)
	{
		FENode& node = mesh.Node(el.m_node[j/3]);
		nj = j%3;

		switch (nj)
		{
		case 0: node.add(dof_X, dx); node.m_rt.x += dx; break;
		case 1: node.add(dof_Y, dx); node.m_rt.y += dx; break;
		case 2: node.add(dof_Z, dx); node.m_rt.z += dx; break;
		}


		fem.Update();

		zero(f1);
		bd.ElementInternalForce(el, f1);

		switch (nj)
		{
		case 0: node.sub(dof_X, dx); node.m_rt.x -= dx; break;
		case 1: node.sub(dof_Y, dx); node.m_rt.y -= dx; break;
		case 2: node.sub(dof_Z, dx); node.m_rt.z -= dx; break;
		}

		fem.Update();

		for (i=0; i<3*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
	}
}
