#include "stdafx.h"
#include "SkylineSolver.h"

//-----------------------------------------------------------------------------
void colsol_factor(int N, double* values, int* pointers);
void colsol_solve (int N, double* values, int* pointers, double* R);

//-----------------------------------------------------------------------------
bool SkylineSolver::PreProcess()
{
	// We don't need to do any preprocessing for this solver
	return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool SkylineSolver::Factor()
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (m_pA);
	assert(pK);

	colsol_factor(pK->Size(), pK->values(), pK->pointers());
	return true;
}

//-----------------------------------------------------------------------------
bool SkylineSolver::Solve(vector<double>& x, vector<double>& R)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (m_pA);

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
void SkylineSolver::Destroy()
{
	// Nothing to destroy
	LinearSolver::Destroy();
}
