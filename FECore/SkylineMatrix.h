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

//=============================================================================
//! Implements a sparse matrix using the skyline storage

//! This class implements a symmetric sparse matrix where only the values
//! below the skyline are stored.

class FECORE_API SkylineMatrix : public SparseMatrix
{
public:
	SkylineMatrix();
	virtual ~SkylineMatrix();

public: // from SparseMatrix

	void Zero() override;

	void Clear() override;

	void Create(SparseMatrixProfile& mp) override;

	void Assemble(const matrix& ke, const std::vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(const matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override;

	void add(int i, int j, double v) override;

	void set(int i, int j, double v) override;

	// NOTE: This is not implemented yet!
	bool check(int i, int j) override;

	double get(int i, int j) override;

	double diag(int i) override;

	double* values() { return m_pd; }
	int* pointers() { return m_ppointers; }

protected:
	void Create(double* pv, int* pp, int N);

protected:
	double*	m_pd;			//!< matrix values
	int*	m_ppointers;	//!< arrays of indices to diagonal elements
};
