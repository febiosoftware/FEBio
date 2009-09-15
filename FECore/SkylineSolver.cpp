#include "stdafx.h"
#include "SkylineSolver.h"

//-----------------------------------------------------------------------------
void colsol_factor(int N, double* values, int* pointers);
void colsol_solve (int N, double* values, int* pointers, double* R);

//-----------------------------------------------------------------------------
bool SkylineSolver::PreProcess(SparseMatrix& K)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n");
		return false;
	}

	// We don't need to do any preprocessing for this solver

	return LinearSolver::PreProcess(K);
}

//-----------------------------------------------------------------------------
bool SkylineSolver::Factor(SparseMatrix& K)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n");
		return false;
	}
#ifdef DEBUG
	int i, n;
	int* pointers;
	double* values;

	n = pK->Size();
	pointers = pK->pointers();
	values = pK->values();
#endif

	colsol_factor(pK->Size(), pK->values(), pK->pointers());

	return true;
}

//-----------------------------------------------------------------------------
bool SkylineSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& R)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n");
		return false;
	}

	// we need to make a copy of R since colsol overwrites the right hand side vector
	// with the solution
	int neq = pK->Size();
	for (int i=0; i<neq; ++i) x[i] = R[i];
	colsol_solve(pK->Size(), pK->values(), pK->pointers(), x);

	return true;
}

//-----------------------------------------------------------------------------
bool SkylineSolver::Solve(SparseMatrix& K, matrix& x, matrix& R)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n");
		return false;
	}

	// we need to make a copy of R since colsol overwrites the right hand side vector
	// with the solution
	int nrhs = x.rows();
	for (int i=0; i<nrhs; ++i)
	{
		int neq = pK->Size();
		for (int j=0; j<neq; ++j) x[i][j] = R[i][j];

		colsol_solve(K.Size(), pK->values(), pK->pointers(), x[i]);
	}

	return true;
}

//-----------------------------------------------------------------------------
void SkylineSolver::Destroy(SparseMatrix& K)
{
	// Nothing to destroy
	LinearSolver::Destroy(K);
}
