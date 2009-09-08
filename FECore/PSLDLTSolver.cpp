#include "stdafx.h"
#include "PSLDLTSolver.h"

//-----------------------------------------------------------------------------
bool PSLDLTSolver::PreProcess(SparseMatrix& K)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// let's make sure the matrix K is of the correct type
	CompactSymmMatrix* pK = dynamic_cast<CompactSymmMatrix*> (&K);
	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n\n");
		return false;
	}

	// Do the preprocessing
	int nonz;
	double ops;
	PSLDLT_Preprocess(0, pK->Size(), pK->pointers(), pK->indices(), &nonz, &ops);

	return LinearSolver::PreProcess(K);

#endif

}

//-----------------------------------------------------------------------------
bool PSLDLTSolver::Factor(SparseMatrix& K)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// let's make sure the matrix K is of the correct type
	CompactSymmMatrix* pK = dynamic_cast<CompactSymmMatrix*> (&K);
	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n\n");
		return false;
	}

#ifdef DEBUG

	int i, n, nnz, *pointers, *indices;
	double* values;

	n = pK->Size();
	nnz = pK->NonZeroes();
	pointers = pK->pointers();
	indices = pK->indices();
	values = pK->values();
	fprintf(stdout, "\nPointers:");
	for (i=0; i<n; i++) fprintf(stdout, "\n%d", pointers[i]);
	fprintf(stdout, "\nIndices, Values:");
	for (i=0; i<nnz; i++) fprintf(stdout, "\n%d, %g", indices[i], values[i]);
#endif

	// Do the factorization
	PSLDLT_Factor(0, pK->Size(), pK->pointers(), pK->indices(), pK->values());
	return true;

#endif

}

//-----------------------------------------------------------------------------
bool PSLDLTSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& R)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// let's make sure the matrix K is of the correct type
	CompactSymmMatrix* pK = dynamic_cast<CompactSymmMatrix*> (&K);
	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n\n");
		return false;
	}

	// Let's roll !!
	PSLDLT_Solve(0, x, R);

	return true;

#endif
}

//-----------------------------------------------------------------------------
bool PSLDLTSolver::Solve(SparseMatrix& K, matrix& x, matrix& b)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// let's make sure the matrix K is of the correct type
	CompactSymmMatrix* pK = dynamic_cast<CompactSymmMatrix*> (&K);
	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n\n");
		return false;
	}


	// Let's roll !!
	int nrhs = x.rows();
	for (int i=0; i<nrhs; ++i) PSLDLT_Solve(0, x[i], b[i]);

	return true;
#endif
}

//-----------------------------------------------------------------------------
void PSLDLTSolver::Destroy()
{
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
#else
	if (m_bvalid) PSLDLT_Destroy(0);
	LinearSolver::Destroy();
#endif
}
