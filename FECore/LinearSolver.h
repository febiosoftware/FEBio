#pragma once

#include "SparseMatrix.h"
#include "vector.h"
#include "matrix.h"

//-----------------------------------------------------------------------------
// matrix types
#define SPARSE_SYMMETRIC	0
#define SPARSE_UNSYMMETRIC	1

//-----------------------------------------------------------------------------
//! base class for the linear solver classes

//! This class defines several virtual functions that need to be overriden
//! in the derived class

class LinearSolver
{
public:
	LinearSolver() { m_bvalid = false; }
	virtual ~LinearSolver() { Destroy(); }

	virtual bool PreProcess(SparseMatrix& K) { m_bvalid = true; return true; }
	virtual bool Factor(SparseMatrix& K) = 0;
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b) = 0;
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b) = 0;
	virtual void Destroy() { m_bvalid = false; };
	virtual void Destroy(SparseMatrix& K) { m_bvalid = false; };

	//! returns a sparse matrix
	virtual SparseMatrix* GetMatrix(int ntype) = 0;

	//! set the number of threads
	static void SetNumThreads(int n) { m_numthreads = (n>0? n : 1); }

protected:
	bool	m_bvalid;	// flag indication wether a valid matrix structure is ready

	static int	m_numthreads;	// nr of threads to create
};

