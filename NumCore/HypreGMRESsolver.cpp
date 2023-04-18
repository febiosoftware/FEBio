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
#include "HypreGMRESsolver.h"
#ifdef HYPRE
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>
#include <_hypre_IJ_mv.h>
#include <HYPRE_krylov.h>


class HypreGMRESsolver::Implementation
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

public:
	// control parameters
	int		m_maxiter;
	double	m_tol;
	int		m_print_level;

public:
	Implementation() : A(0)
	{
		m_print_level = 0;
		m_maxiter = 1000;
		m_tol = 1e-7;
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
		HYPRE_IJMatrixDestroy(ij_A);
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
		HYPRE_IJVectorDestroy(ij_b);
		HYPRE_IJVectorDestroy(ij_x);
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
	void allocPrecond()
	{
		// Now set up the AMG preconditioner and specify any parameters
		HYPRE_BoomerAMGCreate(&precond);
//		HYPRE_BoomerAMGSetPrintLevel(imp->precond, 1); /* print amg solution info */
//		HYPRE_BoomerAMGSetCoarsenType(imp->precond, 6);
		HYPRE_BoomerAMGSetCoarsenType(precond, 10); /* HMIS-coarsening */
		HYPRE_BoomerAMGSetInterpType(precond, 6); /* extended+i interpolation */
		HYPRE_BoomerAMGSetPMaxElmts(precond, 4);
		HYPRE_BoomerAMGSetAggNumLevels(precond, 2);
//		HYPRE_BoomerAMGSetOldDefault(precond);
//		HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
		HYPRE_BoomerAMGSetRelaxType(precond, 3); /* hybrid Gauss-Seidel or SOR, forward solve */
		HYPRE_BoomerAMGSetStrongThreshold(precond, 0.5);
		HYPRE_BoomerAMGSetNumSweeps(precond, 1);
//		HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
	//	HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
	}

	// destroy preconditioner
	void destroyPrecond()
	{
		if (precond) HYPRE_BoomerAMGDestroy(precond);
	}

	// allocate solver
	void allocSolver()
	{
		// Create the solver object
		HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		int    restart = 30;
		HYPRE_FlexGMRESSetKDim(solver, restart);
		HYPRE_FlexGMRESSetMaxIter(solver, m_maxiter); /* max iterations */
		HYPRE_FlexGMRESSetTol(solver, m_tol); /* conv. tolerance */
		//	HYPRE_FlexGMRESSetPrintLevel(imp->solver, 2); /* print solve info */
		HYPRE_FlexGMRESSetLogging(solver, 1); /* needed to get run info later */

		// Set the preconditioner
		HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
			(HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);
	}

	// destroy the solver
	void destroySolver()
	{
		HYPRE_ParCSRFlexGMRESDestroy(solver);
	}

	// calculate the preconditioner
	void doPrecond()
	{
		HYPRE_ParCSRFlexGMRESSetup(solver, par_A, par_b, par_x);
	}

	// solve the linear system
	void doSolve(double* x)
	{
		HYPRE_ParCSRFlexGMRESSolve(solver, par_A, par_b, par_x);

		/* Run info - needed logging turned on */
		int    num_iterations;
		double final_res_norm;
		HYPRE_FlexGMRESGetNumIterations(solver, (HYPRE_Int*)&num_iterations);
		HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
		if (m_print_level != 0)
		{
			printf("\n");
			printf("Iterations = %d\n", num_iterations);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* get the local solution */
		int neq = equations();
		HYPRE_IJVectorGetValues(ij_x, neq, (HYPRE_Int*)&ind[0], &x[0]);
	}
};

BEGIN_FECORE_CLASS(HypreGMRESsolver, LinearSolver)
	ADD_PARAMETER(imp->m_print_level, "print_level");
	ADD_PARAMETER(imp->m_maxiter    , "maxiter"    );
	ADD_PARAMETER(imp->m_tol        , "tol"        );
END_FECORE_CLASS();

HypreGMRESsolver::HypreGMRESsolver(FEModel* fem) : LinearSolver(fem), imp(new HypreGMRESsolver::Implementation)
{

}

HypreGMRESsolver::~HypreGMRESsolver()
{
	delete imp; imp = 0;
}

void HypreGMRESsolver::SetPrintLevel(int n)
{
	imp->m_print_level = n;
}

void HypreGMRESsolver::SetMaxIterations(int n)
{
	imp->m_maxiter = n;
}

void HypreGMRESsolver::SetConvergencTolerance(double tol)
{
	imp->m_tol = tol;
}

SparseMatrix* HypreGMRESsolver::CreateSparseMatrix(Matrix_Type ntype)
{
	if (ntype == Matrix_Type::REAL_UNSYMMETRIC)
		return (imp->A = new CRSSparseMatrix(0));
	else
		return 0;
}

bool HypreGMRESsolver::SetSparseMatrix(SparseMatrix* A)
{
	CRSSparseMatrix* K = dynamic_cast<CRSSparseMatrix*>(A);
	if (K == 0) return false;

	if (K->isRowBased() == false) return false;
	if (K->Offset() != 0) return false;

	imp->A = K;

	return true;
}

//! clean up
void HypreGMRESsolver::Destroy()
{
	// cleanup
	imp->destroyPrecond();
	imp->destroySolver();

	// destroy matrix
	imp->destroyMatrix();

	// Destroy vectors
	imp->destroyVectors();
}

bool HypreGMRESsolver::PreProcess()
{ 
	// make sure data is valid
	if (imp->isValid() == false) return false;

	// create coefficient matrix
	imp->allocMatrix();

	// allocate rhs and solution vectors
	imp->allocVectors();

	return true; 
}

bool HypreGMRESsolver::Factor()
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

bool HypreGMRESsolver::BackSolve(double* x, double* b)
{
	// make sure data is valid
	if (imp->isValid() == false) return false;

	// nr of equations
	int neq = imp->equations();

	// update the vectors
	imp->updateVectors(x, b);

	// solve 
	imp->doSolve(x);

	return true; 
}

#else
BEGIN_FECORE_CLASS(HypreGMRESsolver, LinearSolver)
END_FECORE_CLASS();

HypreGMRESsolver::HypreGMRESsolver(FEModel* fem) : LinearSolver(fem) {}
HypreGMRESsolver::~HypreGMRESsolver() {}
void HypreGMRESsolver::Destroy() {}
void HypreGMRESsolver::SetPrintLevel(int n) {}
void HypreGMRESsolver::SetMaxIterations(int n) {}
void HypreGMRESsolver::SetConvergencTolerance(double tol) {}
bool HypreGMRESsolver::PreProcess() { return false; }
bool HypreGMRESsolver::Factor() { return false; }
bool HypreGMRESsolver::BackSolve(double* x, double* b) { return false; }
SparseMatrix* HypreGMRESsolver::CreateSparseMatrix(Matrix_Type ntype) { return 0; }
bool HypreGMRESsolver::SetSparseMatrix(SparseMatrix* A) { return false; }

#endif // HYPRE
