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
#include "MatrixProfile.h"
#include "MatrixOperator.h"
#include "matrix.h"
#include <vector>

//-----------------------------------------------------------------------------
//! Base class for sparse matrices

//! This is the base class for the sparse matrix classes and defines the interface
//! to the different matrix classes

class FECORE_API SparseMatrix : public MatrixOperator
{
public:
	//! constructor
	SparseMatrix();

	//! destructor
	virtual ~SparseMatrix();

public:
	//! return number of rows
	int Rows() const { return m_nrow; }

	//! return number of columns
	int Columns() const { return m_ncol; }

	//! is the matrix square?
	bool IsSquare() const { return (m_nrow == m_ncol); }

	//! return number of nonzeros
	int NonZeroes() const { return m_nsize; }

public: // functions to be overwritten in derived classes

	//! set all matrix elements to zero
	virtual void Zero() = 0;

	//! Create a sparse matrix from a sparse-matrix profile
	virtual void Create(SparseMatrixProfile& MP) = 0;

	//! assemble a matrix into the sparse matrix
	virtual void Assemble(const matrix& ke, const std::vector<int>& lm) = 0;

	//! assemble a matrix into the sparse matrix
	virtual void Assemble(const matrix& ke, const std::vector<int>& lmi, const std::vector<int>& lmj) = 0;

	//! check if an entry was allocated
	virtual bool check(int i, int j) = 0;

	//! set entry to value
	virtual void set(int i, int j, double v) = 0;

	//! add value to entry
	virtual void add(int i, int j, double v) = 0;

	//! retrieve value
	virtual double get(int i, int j) { return 0; }

	//! get the diagonal value
	virtual double diag(int i) = 0;

	//! release memory for storing data
	virtual void Clear();

	//! scale matrix
	virtual void scale(const std::vector<double>& L, const std::vector<double>& R);

public:
	//! multiply with vector
	bool mult_vector(double* x, double* r) override { assert(false); return false; }

public:
	// NOTE: The following functions are only used by the compact matrices, but I need to be able to override them
	// for the JFNKMatrix so I've moved them here. 
	virtual double* Values() { return 0; }
	virtual int*    Indices() { return 0; }
	virtual int*    Pointers() { return 0; }
	virtual int     Offset() const { return 0; }

protected:
	// NOTE: These values are set by derived classes
	int	m_nrow, m_ncol;		//!< dimension of matrix
	int	m_nsize;			//!< number of nonzeroes (i.e. matrix elements actually allocated)
};
