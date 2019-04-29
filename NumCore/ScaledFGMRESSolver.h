/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FGMRES_ILU0_Solver.h"
#include "CompactUnSymmMatrix.h"

class ScaledFGMRESSolver : public FGMRES_ILU0_Solver
{
public:
	ScaledFGMRESSolver(FEModel* fem) : FGMRES_ILU0_Solver(fem)
	{
		m_k = 1.0;
	}

	bool Factor() override
	{
		CRSSparseMatrix* A = dynamic_cast<CRSSparseMatrix*>(GetSparseMatrix());
		if (A == nullptr) return false;

		int n0 = m_part[0];
		int n1 = m_part[1];

		// we want to multply all columns of partition 2 by 1/k
		for (CRSSparseMatrix::Iterator it(A); it.valid(); it.next())
		{
			MatrixItem item = it.get();
			if (item.col >= n0)
			{
				it.set(item.val / m_k);
			}
		}

		return FGMRES_ILU0_Solver::Factor();
	}

	bool BackSolve(double* x, double* b) override
	{
		bool ret = FGMRES_ILU0_Solver::BackSolve(x, b);
		if (ret == false) return false;

		int n0 = m_part[0];
		int n1 = m_part[1];
		for (size_t i = n0; i < n0+n1; ++i) x[i] /= m_k;

		return true;
	}

	void SetScaleFactor(double k) { m_k = k; }

protected:
	double			m_k;		//!< scale parameter
};
