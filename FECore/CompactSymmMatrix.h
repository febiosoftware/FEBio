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
#include "CompactMatrix.h"
#include "fecore_api.h"

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format (i.e. column major, lower triangular compact).

//! This class also assumes the matrix is symmetric and therefor only stores
//! the lower triangular matrix

class FECORE_API CompactSymmMatrix : public CompactMatrix
{
public:
	//! class constructor
	CompactSymmMatrix(int offset = 0);

	//! Create the matrix structure from the SparseMatrixProfile.
	void Create(SparseMatrixProfile& mp) override;

	//! Assemble an element matrix into the global matrix
	void Assemble(const matrix& ke, const std::vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(const matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override;

	//! add a matrix item
	void add(int i, int j, double v) override;

	//! set matrix item
	void set(int i, int j, double v) override;

	//! get a matrix item
	double get(int i, int j) override;

	// alternative access
	double operator ()(int i, int j) { return get(i, j); }

	//! return the diagonal component
	double diag(int i) override { return m_pd[m_ppointers[i] - m_offset]; }

	//! multiply with vector
	bool mult_vector(double* x, double* r) override;

	//! see if a matrix element is defined
	bool check(int i, int j) override;

	//! is the matrix symmetric or not
	bool isSymmetric() override { return true; }

	//! this is a column based format
	bool isRowBased() override { return false; }

	//! calculate the inf norm
	double infNorm() const override;

	//! calculate the one norm
	double oneNorm() const override;

	//! do row (L) and column (R) scaling
	void scale(const std::vector<double>& L, const std::vector<double>& R) override;
};
