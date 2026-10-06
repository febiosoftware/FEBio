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
#include <FECore/FENLConstraint.h>
#include <FECore/log.h>
#include <FEBioMech/FEMechModel.h>
#include <FEBioMech/FERigidBody.h>
#include <FEBioMech/FESolidSolver2.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

//-----------------------------------------------------------------------------
FEStiffnessDiagnostic::FEStiffnessDiagnostic(FEModel* fem) : FECoreTask(fem)
{
	m_fp = nullptr;
	m_writeMatrix = false;
	m_nmax = -1;
}

//-----------------------------------------------------------------------------
// Initialize the diagnostic. In this function we build the FE model depending
// on the scenario.
bool FEStiffnessDiagnostic::Init(const char* szarg)
{
	if (szarg && szarg[0])
	{
		if (strcmp(szarg, "v") == 0) m_writeMatrix = true;
		else { m_nmax = atoi(szarg); m_writeMatrix = true; }
	}
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

//	fem.AddCallback(stiffness_diagnostic_cb, CB_MATRIX_REFORM, (void*)this);
	fem.AddCallback(stiffness_diagnostic_cb, CB_QUASIN_CONVERGED, (void*)this);

	// create a file name for the log file
	string logfile("diagnostic.log");
	m_fp = fopen(logfile.c_str(), "wt");
	fprintf(m_fp, "FEBio Stiffness Diagnostics:\n");
	fprintf(m_fp, "============================\n");
	fflush(m_fp);

	fem.BlockLog();
	bool bret = fem.Solve();
	fem.UnBlockLog();
	if (bret == false)
	{
		feLogError("FEBio error terminated. Aborting diagnostic.\n");
		fprintf(m_fp, "FEBio error terminated. Diagnostic aborted.\n");
		fclose(m_fp);
		m_fp = nullptr;
		return false;
	}

	fprintf(m_fp, "diagnostic completed.\n");

	fclose(m_fp);
	m_fp = nullptr;

	return true;
}

//-----------------------------------------------------------------------------
// Compares the assembled stiffness matrix with a central finite-difference
// approximation of the derivative of the global residual, evaluated at the
// converged state of each time step. In addition to the overall maximum error,
// a summary is written for each block of the matrix, where blocks are defined
// by the names of the degrees of freedom associated with the rows and columns
// (e.g. x, p, c1, ...). This helps locate inconsistent tangent terms.
bool FEStiffnessDiagnostic::Diagnose()
{
	FEModel* fem = GetFEModel();
	FEMechModel* mech = dynamic_cast<FEMechModel*>(fem);

	FEAnalysis* step = fem->GetCurrentStep();
	if (step == nullptr) return false;

	// any Newton solver will do (e.g. solid, biphasic, multiphasic, ...)
	FENewtonSolver* nlsolve = dynamic_cast<FENewtonSolver*>(step->GetFESolver());
	if (nlsolve == nullptr)
	{
		fprintf(m_fp, "stiffness diagnostic requires a Newton solver.\n");
		return false;
	}

	SparseMatrix* pA = nlsolve->m_pK->GetSparseMatrixPtr();
	if (pA == nullptr) return false;

	// re-evaluate the stiffness matrix at the converged state, since the matrix in
	// memory was evaluated at the start of the last iteration.
	nlsolve->m_pK->Zero();
	std::fill(nlsolve->m_Fd.begin(), nlsolve->m_Fd.end(), 0.0);
	if (nlsolve->StiffnessMatrix() == false)
	{
		fprintf(m_fp, "failed to evaluate the stiffness matrix.\n");
		return false;
	}

	const int neq = pA->Rows();

	// label each equation with the name of its degree of freedom and its node
	DOFS& dofs = fem->GetDOFS();
	std::vector<std::string> eqname(neq, "other");
	std::vector<int> eqnode(neq, -1);

	vector<int> bc(neq, 0);
	int nmax = -1;
	FEMesh& mesh = fem->GetMesh();
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j < (int)node.m_ID.size(); ++j)
		{
			int n = node.m_ID[j];
			if ((n >= 0) && (n < neq))
			{
				const char* sz = dofs.GetDOFName(j);
				eqname[n] = (sz ? sz : "?");
				eqnode[n] = i + 1;
			}
		}

		if (node.m_rid < 0)
		{
			for (int j = 0; j < node.m_ID.size(); ++j)
			{
				int n = node.m_ID[j];
				if (n >= 0) bc[n] = 1;
				if (n > nmax) nmax = n;
			}
		}
		else
		{
			for (int j = 0; j < node.m_ID.size(); ++j)
			{
				int n = -node.m_ID[j]-2;
				if (n >= 0) bc[n] = 1;
				if (n > nmax) nmax = n;
			}
		}
	}

	if (mech)
	{
		for (int i = 0; i < mech->RigidBodies(); ++i)
		{
			FERigidBody& rb = *mech->GetRigidBody(i);
			for (int j = 0; j < 6; ++j)
			{
				int n = rb.m_LM[j];
				if ((n >= 0) && (n < neq)) eqname[n] = "rigid";
				if (n >= 0) bc[n] = 1;
				if (n > nmax) nmax = n;
			}
		}
	}

	if (nmax < neq)
	{
		// these are probably lagrange multiplier dofs
		for (int i = nmax + 1; i < neq; ++i) bc[i] = 1;
	}

	// assign a block index to each dof name
	std::vector<std::string> names;
	std::vector<int> eqblk(neq, 0);
	for (int i = 0; i < neq; ++i)
	{
		int k = -1;
		for (int l = 0; l < (int)names.size(); ++l) if (names[l] == eqname[i]) { k = l; break; }
		if (k < 0) { names.push_back(eqname[i]); k = (int)names.size() - 1; }
		eqblk[i] = k;
	}
	const int nb = (int)names.size();

	// block statistics
	struct BlockStats {
		double kmax = 0;	// max |K_fd|
		double emax = 0;	// max |K - K_fd|
		double e2 = 0;		// sum of squared errors
		double k2 = 0;		// sum of squared K_fd
		int imax = -1, jmax = -1;
	};
	std::vector<BlockStats> blk(nb*nb);

	// current solution (used to scale the perturbations)
	std::vector<double> U(neq, 0.0);
	for (int i = 0; i < neq; ++i)
	{
		if (i < (int)nlsolve->m_Ut.size()) U[i] += nlsolve->m_Ut[i];
		if (i < (int)nlsolve->m_Ui.size()) U[i] += nlsolve->m_Ui[i];
	}

	const double eps = 1e-6;
	int nreq = (m_nmax <= 0 ? neq : m_nmax);
	if (nreq > neq) nreq = neq;

	double max_val = 0, max_err = 0.0;
	int i_max = -1, j_max = -1;
	std::cerr << "\nstarting diagnostic:\nprogress:";
	int pct = 0;

	std::vector<double> u(neq, 0), Rp(neq, 0), Rm(neq, 0);
	for (int j = 0; j < nreq; ++j)
	{
		int new_pct = (100 * j) / nreq;
		if (pct != new_pct) {
			std::cerr << (((new_pct % 10) == 0) ? "+" : "-");
			pct = new_pct;
		}
		if (bc[j] == 0) continue;

		// central difference
		double h = eps*(1.0 + fabs(U[j]));
		std::fill(u.begin(), u.end(), 0.0);
		u[j] = h;
		nlsolve->Update(u);
		std::fill(Rp.begin(), Rp.end(), 0.0);
		nlsolve->Residual(Rp);
		u[j] = -h;
		nlsolve->Update(u);
		std::fill(Rm.begin(), Rm.end(), 0.0);
		nlsolve->Residual(Rm);

		for (int i = 0; i < nreq; ++i)
		{
			if (bc[i] == 0) continue;

			// note that we flip the sign on ka.
			// this is because febio actually calculates the negative of the residual
			double ka_ij = -(Rp[i] - Rm[i]) / (2*h);
			double kt_ij = pA->get(i, j);

			if (fabs(kt_ij) > max_val) max_val = fabs(kt_ij);

			double err = fabs(kt_ij - ka_ij);
			if (err > max_err)
			{
				max_err = err;
				i_max = i;
				j_max = j;
			}

			BlockStats& b = blk[eqblk[i]*nb + eqblk[j]];
			if (fabs(ka_ij) > b.kmax) b.kmax = fabs(ka_ij);
			if (err > b.emax) { b.emax = err; b.imax = i; b.jmax = j; }
			b.e2 += err*err;
			b.k2 += ka_ij*ka_ij;

			if (m_writeMatrix)
			{
				fprintf(m_fp, "%d, %d : %lg, %lg (%lg)\n", i, j, kt_ij, ka_ij, err);
			}
		}
	}
	std::cerr << "\n";

	// let's make sure we leave the model in a consistent state
	std::fill(u.begin(), u.end(), 0.0);
	nlsolve->Update(u);
	nlsolve->Residual(Rp);

	double t = fem->GetCurrentTime();
	fprintf(m_fp, "\n=== time = %lg ===\n", t);
	printf("Max abs. value: %lg\n", max_val);
	fprintf(m_fp, "Max abs. value: %lg\n", max_val);
	if (max_val == 0) max_val = 1;

	printf("Max error: %lg (%d, %d)\n", max_err, i_max, j_max);
	printf("Max rel. error: %lg (%d, %d)\n", max_err / max_val, i_max, j_max);
	fprintf(m_fp, "Max error: %lg (%d, %d)\n", max_err, i_max, j_max);
	fprintf(m_fp, "Max rel. error: %lg (%d, %d)\n", max_err / max_val, i_max, j_max);

	// block summary
	fprintf(m_fp, "\nblock summary (row dof / column dof):\n");
	fprintf(m_fp, "%-8s %-8s %12s %12s %12s   %s\n", "row", "col", "max|Kfd|", "max err", "rel.err(F)", "worst entry: row eq (node), col eq (node)");
	for (int bi = 0; bi < nb; ++bi)
		for (int bj = 0; bj < nb; ++bj)
		{
			BlockStats& b = blk[bi*nb + bj];
			if ((b.kmax == 0) && (b.emax == 0)) continue;
			double rel = (b.k2 > 0 ? sqrt(b.e2 / b.k2) : (b.e2 > 0 ? 1.0 : 0.0));
			fprintf(m_fp, "%-8s %-8s %12.4le %12.4le %12.4le   %d (%d), %d (%d)\n",
				names[bi].c_str(), names[bj].c_str(), b.kmax, b.emax, rel,
				b.imax, (b.imax >= 0 ? eqnode[b.imax] : -1), b.jmax, (b.jmax >= 0 ? eqnode[b.jmax] : -1));
		}
	fflush(m_fp);

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
