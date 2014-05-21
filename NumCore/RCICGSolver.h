// This class implements an interface to the RCI CG iterative solver from the MKL math library.
#pragma once

#include "FECore/LinearSolver.h"
#include "CompactMatrix.h"

class RCICGSolver : public LinearSolver
{
public:
	RCICGSolver();
	virtual bool PreProcess();
	virtual bool Factor();
	virtual bool BackSolve(vector<double>& x, vector<double>& b);
	virtual void Destroy();

	virtual SparseMatrix* CreateSparseMatrix(Matrix_Type ntype);

private:
	CompactSymmMatrix*	m_pA;
};
