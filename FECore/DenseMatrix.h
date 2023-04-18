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

namespace FECore {

//=============================================================================
//! This class implements a full matrix
//! that is a matrix that stores all its elements.

class FECORE_API DenseMatrix : public SparseMatrix
{
public:
	// con/de-structor
	DenseMatrix();
	~DenseMatrix();

	// create a matrix of particular size
	void Create(int rows, int cols);

	// retrieve matrix data
	double& operator () (int i, int j) { return m_pr[i][j]; }

public:
	// zero matrix elements
	void Zero() override;

	// create a matrix from a spares matrix profile
	void Create(SparseMatrixProfile& mp) override;

	// assemble matrix into sparse matrix
	void Assemble(const matrix& ke, const std::vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(const matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override;

	// clear all data
	void Clear() override;

	bool check(int i, int j) override { return true; }
	void add(int i, int j, double v) override { m_pr[i][j] += v; }
	void set(int i, int j, double v) override { m_pr[i][j] = v; }
	double diag(int i) override { return m_pr[i][i]; }

protected:
	double*		m_pd;	//!< matrix values
	double**	m_pr;	//!< pointers to rows
};

}
