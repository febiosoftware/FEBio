#pragma once
#include "FGMRES_ILU0_Solver.h"
#include <FECore/CompactUnSymmMatrix.h>

class ScaledFGMRESSolver : public FGMRES_ILU0_Solver
{
public:
	ScaledFGMRESSolver(FEModel* fem) : FGMRES_ILU0_Solver(fem)
	{
		m_k = 1.0;
	}

	void SetPartitions(const vector<int>& part) override
	{
		m_npart = part;
	}

	bool Factor() override
	{
		CRSSparseMatrix* A = dynamic_cast<CRSSparseMatrix*>(GetSparseMatrix());
		if (A == nullptr) return false;

		int n0 = m_npart[0];
		int n1 = m_npart[1];

		// we want to multply all columns of partition 2 by 1/k
		double* pd = A->Values();
		int* pointers = A->Pointers();
		int* indices = A->Indices();
		int offset = A->Offset();

		// loop over all rows
		int N = n0 + n1;
		for (int i = 0; i<N; ++i)
		{
			double* pv = pd + (pointers[i] - offset);
			int* pi = indices + (pointers[i] - offset);
			int n = pointers[i + 1] - pointers[i];
			for (int j = 0; j < n; ++j) 
				if (pi[j] - offset >= n0)
				{
					pv[j] /= m_k;
				}
		}

		return FGMRES_ILU0_Solver::Factor();
	}

	bool BackSolve(double* x, double* b) override
	{
		bool ret = FGMRES_ILU0_Solver::BackSolve(x, b);
		if (ret == false) return false;

		int n0 = m_npart[0];
		int n1 = m_npart[1];
		for (size_t i = n0; i < n0+n1; ++i) x[i] /= m_k;

		return true;
	}

	void SetScaleFactor(double k) { m_k = k; }

protected:
	vector<int>		m_npart;	//!< where to partition the matrix
	double			m_k;		//!< scale parameter
};
