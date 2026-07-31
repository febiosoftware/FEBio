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
#include "fecore_api.h"

//! This class stores a sparse matrix in Harwell-Boeing format (i.e. column major, lower triangular compact).

//! This class also assumes the matrix is symmetric and therefor only stores
//! the lower triangular matrix

class FECORE_API CompactSymmMatrix64 : public SparseMatrix
{
public:
	//! class constructor
	CompactSymmMatrix64(int offset = 0);
	~CompactSymmMatrix64();

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
	double diag(int i) override { return m_pvalues[m_ppointers[i] - m_offset]; }

	//! multiply with vector

	//! see if a matrix element is defined
	bool check(int i, int j) override;

	void Zero() override;

public:
	double* Values() override { return m_pvalues; }
	long long* Indices64() { return m_pindices; }
	long long* Pointers64() { return m_ppointers; }
	int Offset() const override { return (int)m_offset; }

private:
	int* Indices() override { return nullptr; }
	int* Pointers() override { return nullptr; }
	bool mult_vector(double* x, double* r) override { return false; }

	void Clear();

private:
	double* m_pvalues;			//!< matrix values
	long long* m_pindices;		//!< row indices
	long long* m_ppointers;	//!< pointers
	long long m_offset;		//!< adjust array indices for fortran arrays
};
