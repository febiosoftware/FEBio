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
#include "Hypre_PCG_AMG.h"
#include <FECore/FEModel.h>
#include <FECore/FESolver.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#ifdef HYPRE
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>
#include <_hypre_IJ_mv.h>
#include <HYPRE_krylov.h>


class Hypre_PCG_AMG::Implementation
{
public:
	CRSSparseMatrix*	A;	// global matrix

	vector<int>			ind; // indices array

	// Hypre stuff
	HYPRE_IJMatrix		ij_A;
	HYPRE_ParCSRMatrix	par_A;
	HYPRE_Solver		solver;
	HYPRE_Solver		precond;
	HYPRE_IJVector		ij_b, ij_x;
	HYPRE_ParVector		par_b, par_x;

	FEModel*	m_fem;

	int		m_iters;

public:
	// control parameters
	int		m_maxiter;
	double	m_tol;
	int		m_print_level;
	double	m_amg_tol;

public:
	Implementation() : A(0)
	{
		m_print_level = 0;
		m_maxiter = 1000;
		m_tol = 1e-7;
		m_amg_tol = 1e-7;

		m_iters = 0;

		ij_A = nullptr;
		ij_b = nullptr;
		ij_x = nullptr;
		solver = nullptr;
		precond = nullptr;
		par_A = nullptr;
		par_b = nullptr;
		par_x = nullptr;
	}

	bool isValid() const
	{
		return (A != 0);
	}

	int equations() const { return (A ? A->Rows() : 0); }

	// Allocate stiffness matrix
	void allocMatrix()
	{
		int neq = equations();

		// Create an empty matrix object
		int ret = 0;
		ret = HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, neq - 1, 0, neq - 1, &ij_A);

		// set the matrix object type
		ret = HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
	}

	// destroy stiffness matrix
	void destroyMatrix()
	{
		if (ij_A) HYPRE_IJMatrixDestroy(ij_A);
		ij_A = nullptr;
	}

	// Allocate vectors for rhs and solution
	void allocVectors()
	{
		int neq = equations();
		ind.resize(neq, 0);
		for (int i = 0; i<neq; ++i) ind[i] = i;

		// create the vector object for the rhs
		HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, neq - 1, &ij_b);
		HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);

		// create the vector object for the solution
		HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, neq - 1, &ij_x);
		HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
	}

	// destroy vectors
	void destroyVectors()
	{
		if (ij_b) HYPRE_IJVectorDestroy(ij_b); ij_b = nullptr;
		if (ij_x) HYPRE_IJVectorDestroy(ij_x); ij_x = nullptr;
	}

	// update vectors 
	void updateVectors(double* x, double* b)
	{
		// initialize vectors for changing coefficient values
		HYPRE_IJVectorInitialize(ij_b);
		HYPRE_IJVectorInitialize(ij_x);

		// set the values
		int neq = equations();
		HYPRE_IJVectorSetValues(ij_b, neq, (HYPRE_Int*)&ind[0], b);
		HYPRE_IJVectorSetValues(ij_x, neq, (HYPRE_Int*)&ind[0], x);

		// finialize assembly
		HYPRE_IJVectorAssemble(ij_b);
		HYPRE_IJVectorAssemble(ij_x);

		HYPRE_IJVectorGetObject(ij_b, (void**)&par_b);
		HYPRE_IJVectorGetObject(ij_x, (void**)&par_x);
	}

	// update coefficient matrix
	void updateMatrix()
	{
		int neq = equations();

		// call initialize, after which we can set the matrix coefficients
		HYPRE_Int ret = HYPRE_IJMatrixInitialize(ij_A);

		// set the matrix coefficients
		double* values = A->Values();
		int* indices = A->Indices();
		int* pointers = A->Pointers();
		for (int i = 0; i<neq; ++i)
		{
			const int* cols = indices + pointers[i];
			int ncols = pointers[i + 1] - pointers[i];
			double* vals = values + pointers[i];
			HYPRE_Int nrows = 1;
			ret = HYPRE_IJMatrixSetValues(ij_A, nrows, (HYPRE_Int*)&ncols, (HYPRE_Int*)&i, (HYPRE_Int*)cols, vals);
		}

		// Finalize the matrix assembly
		ret = HYPRE_IJMatrixAssemble(ij_A);

		// get the matrix object for later use
		ret = HYPRE_IJMatrixGetObject(ij_A, (void**)&par_A);
	}

	// allocate preconditioner
	bool allocPrecond()
	{
		// Now set up the AMG preconditioner and specify any parameters
		HYPRE_BoomerAMGCreate(&precond);
//		HYPRE_BoomerAMGSetPrintLevel(precond, 2);			// print solve info
		HYPRE_BoomerAMGSetCoarsenType(precond, 10);			// HMIS-coarsening
		HYPRE_BoomerAMGSetStrongThreshold(precond, 0.5);	// Threshold for choosing weak/strong connections
		HYPRE_BoomerAMGSetMaxRowSum(precond, 1.0);			// Disable dependency weakening based on maximum row sum.
		HYPRE_BoomerAMGSetCycleRelaxType(precond, 13, 1);	// Pre-smoother: forward L1-Gauss-Seidel
		HYPRE_BoomerAMGSetCycleRelaxType(precond, 14, 2);	// Post-smoother: backward L1-Gauss-Seidel
		HYPRE_BoomerAMGSetCycleRelaxType(precond, 9, 3);	// Coarsest grid solver: Gauss elimination
		HYPRE_BoomerAMGSetAggNumLevels(precond, 1);			// One level of aggressive coarsening
		HYPRE_BoomerAMGSetNumPaths(precond, 1);				// Number of paths of length 2 for aggressive coarsening
		HYPRE_BoomerAMGSetInterpType(precond, 6);			// Extended+i interpolation
		HYPRE_BoomerAMGSetMaxIter(precond, m_maxiter);      // Set maximum iterations
		HYPRE_BoomerAMGSetTol(precond, m_amg_tol);			// conv. tolerance

		FESolver* fesolve = m_fem->GetCurrentStep()->GetFESolver();

		// get the dof map
		vector<int> dofMap;
		int nfunc = fesolve->GetActiveDofMap(dofMap);
		if (nfunc == -1) return false;

		// allocate dof map
		// (We need to copy it here since Hypre will deallocate it)
		int neq = (int)dofMap.size();
		int* dof_func = (int*)malloc(neq * sizeof(int));
		for (size_t i = 0; i < neq; ++i) dof_func[i] = dofMap[i];

		printf("\tNumber of functions : %d\n", nfunc);

		// assign to BoomerAMG
		HYPRE_BoomerAMGSetNumFunctions(precond, nfunc);

		// set the dof map
		HYPRE_BoomerAMGSetDofFunc(precond, (HYPRE_Int*)dof_func);

		return true;
	}

	// destroy preconditioner
	void destroyPrecond()
	{
		if (precond) HYPRE_BoomerAMGDestroy(precond);
		precond = nullptr;
	}

	// allocate solver
	void allocSolver()
	{
		// Create the solver object
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetTwoNorm(solver, 1);
		HYPRE_PCGSetTol(solver, m_tol);

		// Set the preconditioner
		HYPRE_ParCSRPCGSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve,
			(HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup, precond);
	}

	// destroy the solver
	void destroySolver()
	{
		if (solver) HYPRE_ParCSRPCGDestroy(solver);
		solver = nullptr;
	}

	// calculate the preconditioner
	void doPrecond()
	{
		HYPRE_ParCSRPCGSetup(solver, par_A, par_b, par_x);
	}

	// solve the linear system
	void doSolve(double* x)
	{
		HYPRE_ParCSRPCGSolve(solver, par_A, par_b, par_x);

		/* Run info - needed logging turned on */
		double final_res_norm;
		HYPRE_ParCSRPCGGetNumIterations(solver, (HYPRE_Int*)&m_iters);
		HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (m_print_level != 0)
		{
			feLogEx(m_fem, "\n");
			feLogEx(m_fem, "Iterations = %d\n", m_iters);
			feLogEx(m_fem, "Final Relative Residual Norm = %e\n", final_res_norm);
			feLogEx(m_fem, "\n");
		}

		/* get the local solution */
		int neq = equations();
		HYPRE_IJVectorGetValues(ij_x, neq, (HYPRE_Int*)&ind[0], &x[0]);
	}
};

BEGIN_FECORE_CLASS(Hypre_PCG_AMG, LinearSolver)
	ADD_PARAMETER(imp->m_print_level, "print_level");
	ADD_PARAMETER(imp->m_maxiter    , "maxiter"    );
	ADD_PARAMETER(imp->m_tol        , "tol"        );
	ADD_PARAMETER(imp->m_amg_tol    , "amg_tol"    );
END_FECORE_CLASS();

Hypre_PCG_AMG::Hypre_PCG_AMG(FEModel* fem) : LinearSolver(fem), imp(new Hypre_PCG_AMG::Implementation)
{
	imp->m_fem = fem;
}

Hypre_PCG_AMG::~Hypre_PCG_AMG()
{
	delete imp; imp = 0;
}

void Hypre_PCG_AMG::SetPrintLevel(int n)
{
	imp->m_print_level = n;
}

void Hypre_PCG_AMG::SetMaxIterations(int n)
{
	imp->m_maxiter = n;
}

void Hypre_PCG_AMG::SetConvergencTolerance(double tol)
{
	imp->m_tol = tol;
}

SparseMatrix* Hypre_PCG_AMG::CreateSparseMatrix(Matrix_Type ntype)
{
	if (ntype == Matrix_Type::REAL_UNSYMMETRIC)
		return (imp->A = new CRSSparseMatrix(0));
	else
		return 0;
}

bool Hypre_PCG_AMG::SetSparseMatrix(SparseMatrix* A)
{
	CRSSparseMatrix* K = dynamic_cast<CRSSparseMatrix*>(A);
	if (K == 0) return false;

	if (K->isRowBased() == false) return false;
	if (K->Offset() != 0) return false;

	imp->A = K;

	return true;
}

//! clean up
void Hypre_PCG_AMG::Destroy()
{
	// cleanup
	imp->destroyPrecond();
	imp->destroySolver();

	// destroy matrix
	imp->destroyMatrix();

	// Destroy vectors
	imp->destroyVectors();
}

bool Hypre_PCG_AMG::PreProcess()
{ 
	// make sure data is valid
	if (imp->isValid() == false) return false;

	// create coefficient matrix
	imp->allocMatrix();

	// allocate rhs and solution vectors
	imp->allocVectors();

	return true; 
}

bool Hypre_PCG_AMG::Factor()
{ 
	// make sure data is valid
	if (imp->isValid() == false) return false;

	// copy matrix values
	imp->updateMatrix();

	// initialize vectors
	int neq = imp->equations();
	vector<double> zero(neq, 0.0);
	imp->updateVectors(&zero[0], &zero[0]);

	// allocate preconditioner (always call before creating solver!)
	imp->allocPrecond();

	// allocate solver
	imp->allocSolver();

	// apply preconditioner
	imp->doPrecond();

	return true;
}

bool Hypre_PCG_AMG::BackSolve(double* x, double* b)
{
	// make sure data is valid
	if (imp->isValid() == false) return false;

	// nr of equations
	int neq = imp->equations();

	// update the vectors
	imp->updateVectors(x, b);

	// solve 
	imp->doSolve(x);

	// update stats
	UpdateStats(imp->m_iters);

	return true; 
}

#else
BEGIN_FECORE_CLASS(Hypre_PCG_AMG, LinearSolver)
END_FECORE_CLASS();

Hypre_PCG_AMG::Hypre_PCG_AMG(FEModel* fem) : LinearSolver(fem) {}
Hypre_PCG_AMG::~Hypre_PCG_AMG() {}
void Hypre_PCG_AMG::Destroy() {}
void Hypre_PCG_AMG::SetPrintLevel(int n) {}
void Hypre_PCG_AMG::SetMaxIterations(int n) {}
void Hypre_PCG_AMG::SetConvergencTolerance(double tol) {}
bool Hypre_PCG_AMG::PreProcess() { return false; }
bool Hypre_PCG_AMG::Factor() { return false; }
bool Hypre_PCG_AMG::BackSolve(double* x, double* b) { return false; }
SparseMatrix* Hypre_PCG_AMG::CreateSparseMatrix(Matrix_Type ntype) { return 0; }
bool Hypre_PCG_AMG::SetSparseMatrix(SparseMatrix* A) { return false; }

#endif // HYPRE
