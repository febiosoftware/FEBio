/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "SparseMatrix.h"

class FENewtonSolver;

//-----------------------------------------------------------------------------
// This is a class that just mimics a sparse matrix.
// It is only used by the JFNK strategy. 
// The only function it implements is the mult_vector.
class FECORE_API JFNKMatrix : public SparseMatrix
{
public:
	enum MultiplyPolicy {
		ZERO_FREE_DOFS,
		ZERO_PRESCRIBED_DOFS
	};

public:
	JFNKMatrix(FENewtonSolver* pns, SparseMatrix* K = 0);

	//! override multiply with vector (Does not use sparse matrix m_K)
	bool mult_vector(double* x, double* r) override;

	//! set the reference residual
	void SetReferenceResidual(std::vector<double>& R0);

	//! set matrix policy
	void SetPolicy(MultiplyPolicy p);

	//! set the forward difference epsilon
	void SetEpsilon(double eps);

public: // these functions use the actual sparse matrix m_K

	//! set all matrix elements to zero
	void Zero() override { m_K->Zero(); }

	//! Create a sparse matrix from a sparse-matrix profile
	void Create(SparseMatrixProfile& MP) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(const matrix& ke, const std::vector<int>& lm) override { m_K->Assemble(ke, lm); }

	//! assemble a matrix into the sparse matrix
	void Assemble(const matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override { m_K->Assemble(ke, lmi, lmj); }

	//! check if an entry was allocated
	bool check(int i, int j) override { return m_K->check(i, j); }

	//! set entry to value
	void set(int i, int j, double v) override { m_K->set(i, j, v); }

	//! add value to entry
	void add(int i, int j, double v) override { m_K->add(i, j, v); }

	//! retrieve value
	double get(int i, int j) override { return m_K->get(i, j); }

	//! get the diagonal value
	double diag(int i) override { return m_K->diag(i); }

	//! release memory for storing data
	void Clear() override { m_K->Clear(); }

	// interface to compact matrices
	double* Values() override { return m_K->Values(); }
	int*    Indices() override { return m_K->Indices(); }
	int*    Pointers() override { return m_K->Pointers(); }
	int     Offset() const override { return m_K->Offset(); }

private:
	bool			m_bauto_eps;	// calculate epsilon automatically
	double			m_eps;		// forward difference epsilon
	SparseMatrix*	m_K;		// the actual sparse matrix (This is only used as a preconditioner and can be null)
	FENewtonSolver*	m_pns;
	std::vector<double>	m_v, m_R;

	std::vector<double>	m_R0;

	std::vector<int>		m_freeDofs, m_prescribedDofs;
	MultiplyPolicy	m_policy;
};
