#pragma once

#include "FECore/LinearSolver.h"
#include "CompactMatrix.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! Implements a linear solver that uses a compact column storage format.

//! This solver can only be used on systems where it is available,
//! such as SGI IRIX, SGI ALTIX, ...

#ifdef PSLDLT
extern "C" {
	void PSLDLT_Ordering(int token, int method);
	void PSLDLT_Preprocess(int, int, int*, int*, int*, double*);
	void PSLDLT_Factor(int, int, int*, int*, double*);
	void PSLDLT_Solve(int, double*, double*);
	void PSLDLT_Destroy(int token);
}
#endif // PSLDLT

class PSLDLTSolver : public LinearSolver
{
public:
	bool PreProcess();
	bool Factor();
	bool BackSolve(vector<double>& x, vector<double>& b);
	void Destroy();

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) { return (m_pA = (ntype == SPARSE_SYMMETRIC? new CompactSymmMatrix() : 0)); }
};
