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
#include <FECore/SparseMatrix.h>
#include <FECore/LinearSolver.h>
#include <FECore/CompactSymmMatrix.h>
#include <FECore/CompactUnSymmMatrix.h>

//-----------------------------------------------------------------------------
// This class implements a diagonally symmetric block-structured matrix. That is
// A matrix for which the diagonal blocks are symmetric, but the off-diagonal
// matrices can be unsymmetric.
class BlockMatrix : public SparseMatrix
{
public:
	struct BLOCK
	{
		int		nstart_row, nend_row;
		int		nstart_col, nend_col;
		CompactMatrix*	pA;

		int Rows() { return nend_row - nstart_row + 1; }
		int Cols() { return nend_col - nstart_col + 1; }

		bool vmult(vector<double>& x, vector<double>& y)
		{
			if (pA) return pA->mult_vector(&x[0], &y[0]);
			else for (size_t i = 0; i < x.size(); ++i) y[i] = 0.0;
			return true;
		}
	};

public:
	BlockMatrix();
	~BlockMatrix();

public:
	//! Partition the matrix into blocks
	void Partition(const vector<int>& part, Matrix_Type mtype, int offset = 1);

public:
	//! Create a sparse matrix from a sparse-matrix profile
	void Create(SparseMatrixProfile& MP) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(const matrix& ke, const std::vector<int>& lm) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(const matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) override;

	//! check if a matrix entry was allocated
	bool check(int i, int j) override;

	//! set entry to value
	void set(int i, int j, double v) override;

	//! add value to entry
	void add(int i, int j, double v) override;

	//! retrieve value
	double get(int i, int j) override;

	//! get the diagonal value
	double diag(int i) override;

	//! release memory for storing data
	void Clear() override;

	//! zero matrix elements
	void Zero() override;

	//! multiply with vector
	bool mult_vector(double* x, double* r) override;

	//! row and column scale
	void scale(const vector<double>& L, const vector<double>& R) override;

public:
	//! return number of blocks
	int Blocks() const { return (int) m_Block.size(); }

	//! get a block
	BLOCK& Block(int i, int j);

	//! find the partition index of an equation number i
	int find_partition(int i);

	//! return number of partitions
	int Partitions() const { return (int) m_part.size() - 1; }

	//! Start equation index of partition i
	int StartEquationIndex(int i) { return m_part[i]; }

	//! number of equations in partition i
	int PartitionEquations(int i) { return m_part[i+1]-m_part[i]; }

protected:
	vector<int>		m_part;		//!< partition list
	vector<BLOCK>	m_Block;	//!< block matrices
};
