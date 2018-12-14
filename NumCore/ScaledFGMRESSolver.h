#pragma once
#include "FGMRES_ILU0_Solver.h"
#include <FECore/CompactUnSymmMatrix.h>

class ScaledFGMRESSolver : public FGMRESSolver
{
public:
	ScaledFGMRESSolver(FEModel* fem) : FGMRESSolver(fem)
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
/*
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
			for (int j = 0; j < n; ++j) if (pi[j] - offset >= n0) pv[j] /= m_k;
		}
*/
		return FGMRESSolver::Factor();
	}

	void mult_vector(double* x, double* y) override
	{
		int n0 = m_npart[0];
		int n1 = m_npart[1];
		int N = n0 + n1;
		vector<double> tmp(N);
		for (size_t i = 0; i < n0; ++i) tmp[i] = x[i];
		for (size_t i = n0; i < N; ++i) tmp[i] = x[i] / m_k;
		FGMRESSolver::mult_vector(&tmp[0], y);
	}

	bool BackSolve(double* x, double* b) override
	{
		bool ret = FGMRESSolver::BackSolve(x, b);
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
