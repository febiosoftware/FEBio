#include "stdafx.h"
#include "HypreGMRESsolver.h"
#ifdef HYPRE
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>


class HypreGMRESsolver::Implementation
{
public:
	CompactUnSymmMatrix*	A;	// global matrix

	// Hypre stuff
	HYPRE_IJMatrix		ij_A;
	HYPRE_ParCSRMatrix	par_A;

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

	int equations() const { return (A ? A->Size() : 0); }
};

HypreGMRESsolver::HypreGMRESsolver() : imp(new HypreGMRESsolver::Implementation)
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
		return (imp->A = new CompactUnSymmMatrix(0, true));
	else
		return 0;
}

bool HypreGMRESsolver::PreProcess()
{ 
	// make sure data is valid
	if (imp->isValid() == false) return false;

	// get the size of the matrix (i.e. nr of rows, cols)
	int neq = imp->equations();

	// Create an empty matrix object
	int ret = 0;
	ret = HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, neq - 1, 0, neq - 1, &imp->ij_A);

	// set the matrix object type
	ret = HYPRE_IJMatrixSetObjectType(imp->ij_A, HYPRE_PARCSR);

	return true; 
}

bool HypreGMRESsolver::Factor()
{ 
	// make sure data is valid
	if (imp->isValid() == false) return false;

	// get the size of the matrix (i.e. nr of rows, cols)
	int neq = imp->equations();

	// call initialize, after which we can set the matrix coefficients
	int ret = HYPRE_IJMatrixInitialize(imp->ij_A);

	// set the matrix coefficients
	double* values = imp->A->Values();
	int* indices = imp->A->Indices();
	int* pointers = imp->A->Pointers();
	for (int i = 0; i<neq; ++i)
	{
		const int* cols = indices + pointers[i];
		int ncols = pointers[i + 1] - pointers[i];
		double* vals = values + pointers[i];
		ret = HYPRE_IJMatrixSetValues(imp->ij_A, 1, &ncols, &i, cols, vals);
	}

	// Finalize the matrix assembly
	ret = HYPRE_IJMatrixAssemble(imp->ij_A);

	// get the matrix object for later use
	ret = HYPRE_IJMatrixGetObject(imp->ij_A, (void**)&imp->par_A);

	return true;
}

bool HypreGMRESsolver::BackSolve(vector<double>& x, vector<double>& b)
{
	// make sure data is valid
	if (imp->isValid() == false) return false;

	// nr of equations
	int neq = imp->equations();

	// create the vector object for the rhs
	HYPRE_IJVector	ij_b;
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, neq - 1, &ij_b);
	HYPRE_IJVectorSetObjectType(ij_b, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(ij_b);

	// create the vector object for the solution
	HYPRE_IJVector	ij_x;
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, neq - 1, &ij_x);
	HYPRE_IJVectorSetObjectType(ij_x, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(ij_x);

	// set the values
	vector<int> indices(neq);
	for (int i=0; i<neq; ++i) indices[i] = i;
	HYPRE_IJVectorSetValues(ij_b, neq, &indices[0], &b[0]);
	HYPRE_IJVectorSetValues(ij_x, neq, &indices[0], &x[0]);

	// finialize assembly
	HYPRE_IJVectorAssemble(ij_x);
	HYPRE_IJVectorAssemble(ij_b);

	// retrieve the vector objects
	HYPRE_ParVector par_b;
	HYPRE_IJVectorGetObject(ij_b, (void**)&par_b);

	HYPRE_ParVector par_x;
	HYPRE_IJVectorGetObject(ij_x, (void**)&par_x);

	int    num_iterations;
	double final_res_norm;
	int    restart = 30;
	int    modify = 1;


	// Create the solver object
	HYPRE_Solver solver;
	HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_FlexGMRESSetKDim(solver, restart);
	HYPRE_FlexGMRESSetMaxIter(solver, imp->m_maxiter); /* max iterations */
	HYPRE_FlexGMRESSetTol(solver, imp->m_tol); /* conv. tolerance */
//	HYPRE_FlexGMRESSetPrintLevel(solver, 2); /* print solve info */
	HYPRE_FlexGMRESSetLogging(solver, 1); /* needed to get run info later */


	/* Now set up the AMG preconditioner and specify any parameters */
	HYPRE_Solver precond;
	HYPRE_BoomerAMGCreate(&precond);
//	HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
	HYPRE_BoomerAMGSetCoarsenType(precond, 6);
	HYPRE_BoomerAMGSetOldDefault(precond);
	HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
	HYPRE_BoomerAMGSetNumSweeps(precond, 1);
	HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
	HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

	/* Set the FlexGMRES preconditioner */
	HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
		(HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);


	/* Now setup and solve! */
	HYPRE_ParCSRFlexGMRESSetup(solver, imp->par_A, par_b, par_x);
	HYPRE_ParCSRFlexGMRESSolve(solver, imp->par_A, par_b, par_x);

	/* Run info - needed logging turned on */
	HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
	if (imp->m_print_level != 0)
	{
		printf("\n");
		printf("Iterations = %d\n", num_iterations);
		printf("Final Relative Residual Norm = %e\n", final_res_norm);
		printf("\n");
	}

	/* Destory solver and preconditioner */
	HYPRE_ParCSRFlexGMRESDestroy(solver);
	HYPRE_BoomerAMGDestroy(precond);

	/* get the local solution */
	HYPRE_IJVectorGetValues(ij_x, neq, &indices[0], &x[0]);

	return true; 
}

#else
HypreGMRESsolver::HypreGMRESsolver(){}
HypreGMRESsolver::~HypreGMRESsolver() {}
void HypreGMRESsolver::SetPrintLevel(int n) {}
void HypreGMRESsolver::SetMaxIterations(int n) {}
void HypreGMRESsolver::SetConvergencTolerance(double tol) {}
bool HypreGMRESsolver::PreProcess() { return false; }
bool HypreGMRESsolver::Factor() { return false; }
bool HypreGMRESsolver::BackSolve(vector<double>& x, vector<double>& b) { return false; }
SparseMatrix* HypreGMRESsolver::CreateSparseMatrix(Matrix_Type ntype) { return 0; }

#endif // HYPRE
