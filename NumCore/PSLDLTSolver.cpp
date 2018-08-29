#include "stdafx.h"
#include "PSLDLTSolver.h"

//-----------------------------------------------------------------------------
PSLDLTSolver::PSLDLTSolver() : m_pA(0)
{
}

//-----------------------------------------------------------------------------
SparseMatrix* PSLDLTSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	return (m_pA = (ntype == REAL_SYMMETRIC? new CompactSymmMatrix() : 0)); 
}

//-----------------------------------------------------------------------------
bool PSLDLTSolver::PreProcess()
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// Do the preprocessing
	int nonz;
	double ops;
	PSLDLT_Preprocess(0, m_pA->Size(), m_pA->pointers(), m_pA->indices(), &nonz, &ops);

	return LinearSolver::PreProcess();
#endif
}


//-----------------------------------------------------------------------------
bool PSLDLTSolver::Factor()
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

#ifdef DEBUG

	int i, n, nnz, *pointers, *indices;
	double* values;

	n = m_pA->Size();
	nnz = m_pA->NonZeroes();
	pointers = m_pA->pointers();
	indices = m_pA->indices();
	values = m_pA->values();
	fprintf(stdout, "\nPointers:");
	for (i=0; i<n; i++) fprintf(stdout, "\n%d", pointers[i]);
	fprintf(stdout, "\nIndices, Values:");
	for (i=0; i<nnz; i++) fprintf(stdout, "\n%d, %g", indices[i], values[i]);
#endif

	// Do the factorization
	PSLDLT_Factor(0, m_pA->Size(), m_pA->pointers(), m_pA->indices(), m_pA->values());
	return true;
#endif
}

//-----------------------------------------------------------------------------
bool PSLDLTSolver::BackSolve(vector<double>& x, vector<double>& R)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// Let's roll !!
	PSLDLT_Solve(0, &x[0], &R[0]);

	return true;

#endif
}

//-----------------------------------------------------------------------------
void PSLDLTSolver::Destroy()
{
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
#else
	LinearSolver::Destroy();
#endif
}
