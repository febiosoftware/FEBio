#pragma once
#include <FECore/LinearSolver.h>

class BoomerAMGSolver : public LinearSolver
{
	class Implementation;

public:
	BoomerAMGSolver(FEModel* fem);
	~BoomerAMGSolver();

public:
	void SetPrintLevel(int printLevel);
	void SetMaxIterations(int maxIter);
	void SetConvergenceTolerance(double tol);
	void SetMaxLevels(int levels);
	void SetCoarsenType(int coarsenType);
	void SetNumFunctions(int funcs);

public:
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
	bool SetSparseMatrix(SparseMatrix* pA) override;

	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* y) override;
	void Destroy() override;

private:
	Implementation*	imp;
};
