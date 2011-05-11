//! This class implements the Watson Sparse Matrix Package.

//! The WSMP solver requires a license file.
//! Documentation can be found at:
//!	http://www-users.cs.umn.edu/~agupta/wsmp

#pragma once

#include "FECore/SparseMatrix.h"
#include "FECore/LinearSolver.h"
#include "CompactMatrix.h"
#include "FECore/vector.h"
#include "FECore/matrix.h"

using namespace FECore;


	/* WSMP Fortran prototypes */
#ifdef WSMP
extern "C"
{
	void wsetmaxthrds_(int *);

	void wsmp_initialize_();

	void wssmp_(int *, int *, int *, double *, double *, int *, int *, double *,
		int *, int *, double *, int *, int *, int *, double *);

	void wsmp_clear_();
}
#endif //WSMP

class WSMPSolver : public LinearSolver
{
public:
	bool PreProcess();
	bool Factor();
	bool BackSolve(vector<double>& x, vector<double>& b);
	void Destroy();

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) { return (m_pA = (ntype == SPARSE_SYMMETRIC? new CompactSymmMatrix(1) : 0)); }

protected:
	/* WSMP control parameters */
	int m_iparm[64];
	double m_dparm[64];

	/* Matrix data */
	int m_n, m_nnz;
	vector<int> m_perm, m_invp;
	vector<double> m_b;
};
