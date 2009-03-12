//! This class implements the Watson Sparse Matrix Package.

//! The WSMP solver requires a license file.
//! Documentation can be found at:
//!	http://www-users.cs.umn.edu/~agupta/wsmp

#pragma once

#include "SparseMatrix.h"
#include "LinearSolver.h"
#include "vector.h"
#include "matrix.h"

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
	virtual bool PreProcess(SparseMatrix& K);
	virtual bool Factor(SparseMatrix& K);
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	virtual void Destroy();

	virtual SparseMatrix* GetMatrix() { return new CompactMatrix(1); }

protected:
	/* WSMP control parameters */
	int m_iparm[64];
	double m_dparm[64];

	/* Matrix data */
	int m_n, m_nnz;
	vector<int> m_perm, m_invp;
	vector<double> m_b;
};

