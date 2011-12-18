#pragma once

#include "SparseMatrix.h"
#include "NumCore/matrix.h"
#include <vector>
using namespace std;

namespace FECore {

//-----------------------------------------------------------------------------
// matrix types
enum Matrix_Type {
	SPARSE_SYMMETRIC,
	SPARSE_UNSYMMETRIC
};

//-----------------------------------------------------------------------------
//! base class for the linear solver classes

//! This class defines several virtual functions that need to be overriden
//! in the derived class

class LinearSolver
{
public:
	LinearSolver() { m_bvalid = false; m_pA = 0; }
	virtual ~LinearSolver() { Destroy(); }

	virtual bool PreProcess() { m_bvalid = true; return true; }
	virtual bool Factor() = 0;
	virtual bool BackSolve(vector<double>& x, vector<double>& b) = 0;
	virtual void Destroy() { m_bvalid = false; };

	//! returns a pointer to the sparse matrix
	SparseMatrix* GetMatrix() { return m_pA; };

	//! set the number of threads
	static void SetNumThreads(int n) { m_numthreads = (n>0? n : 1); }

	// create the sparse matrix
	virtual SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) = 0;

protected:
	bool	m_bvalid;	// flag indication wether a valid matrix structure is ready

	SparseMatrix*	m_pA;	// the matrix that stores the coefficients

	static int	m_numthreads;	// nr of threads to create
};

} // namespace FECore
