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
#include "FEStiffnessDiagnostic.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FENewtonSolver.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
FEStiffnessDiagnostic::FEStiffnessDiagnostic(FEModel* fem) : FECoreTask(fem)
{

}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEStiffnessDiagnostic::Init(const char* szfile)
{
	return GetFEModel()->Init();
}

//-----------------------------------------------------------------------------
bool stiffness_diagnostic_cb(FEModel* fem, unsigned int when, void* pd)
{
	FEStiffnessDiagnostic* diagnostic = (FEStiffnessDiagnostic*)pd;
	return diagnostic->Diagnose();
}

//-----------------------------------------------------------------------------
// Run the tangent diagnostic. After we run the FE model, we calculate 
// the element stiffness matrix and compare that to a finite difference
// of the element residual.
bool FEStiffnessDiagnostic::Run()
{
	// solve the problem
	FEModel& fem = *GetFEModel();

	fem.AddCallback(stiffness_diagnostic_cb, CB_MATRIX_REFORM, (void*)this);

	fem.BlockLog();
	bool bret = fem.Solve();
	fem.UnBlockLog();
	if (bret == false)
	{
		feLogError("FEBio error terminated. Aborting diagnostic.\n");
		return false;
	}

	// create a file name for the log file
	string logfile("diagnostic.log");
	FILE* fp = fopen(logfile.c_str(), "wt");

	fprintf(fp, "FEBio Stiffness Diagnostics Results:\n");
	fprintf(fp, "==================================\n");
	fprintf(fp, "All good!\n");

	fclose(fp);

	return true;
}

//-----------------------------------------------------------------------------
bool FEStiffnessDiagnostic::Diagnose()
{
	FEModel* fem = GetFEModel();

	FEAnalysis* step = fem->GetCurrentStep();
	if (step == nullptr) return false;

	FENewtonSolver* nlsolve = dynamic_cast<FENewtonSolver*>(step->GetFESolver());
	if (nlsolve == nullptr) return false;

	SparseMatrix* pA = nlsolve->m_pK->GetSparseMatrixPtr();
	if (pA == nullptr) return false;

	const double eps = 1e-6;
	int neq = pA->Rows();
	std::vector<double> R0(neq, 0);
	nlsolve->Residual(R0);
	double max_err = 0.0;
	for (int j = 0; j < neq; ++j)
	{
		std::vector<double> u(neq, 0);
		std::vector<double> R(neq, 0);
		u[j] = eps;
		nlsolve->Update(u);
		nlsolve->Residual(R);

		for (int i = 0; i < neq; ++i)
		{
			// note that we flip the sign on ka.
			// this is because febio actually calculates the negative of the residual
			double ka_ij = -(R[i] - R0[i]) / eps;
			double kt_ij = pA->get(i, j);

			double err = fabs(kt_ij - ka_ij);
			if (err > max_err) max_err = err;
		}
	}

	printf("Max error: %lg\n", max_err);

	return true;
}

//-----------------------------------------------------------------------------
// Calculate a finite difference approximation of the derivative of the
// element residual.
void FEStiffnessDiagnostic::deriv_residual(matrix& ke)
{
/*	// get the solver
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
	for (j = 0; j < 3 * N; ++j)
	{
		FENode& node = mesh.Node(el.m_node[j / 3]);
		nj = j % 3;

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

		for (i = 0; i < 3 * N; ++i) ke[i][j] = -(f1[i] - f0[i]) / dx;
	}
*/
}
