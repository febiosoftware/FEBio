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
#include <vector>
#include "fecore_api.h"

// This class represents a sparse matrix in the row-compressed format (3-array format)
class FECORE_API CSRMatrix
{
public:
	// default constructor
	CSRMatrix();

	// create a matrix of given size
	CSRMatrix(int rows, int cols, int noffset = 0);

	// copy constructor
	CSRMatrix(const CSRMatrix& A);

	// Create matrix
	void create(int nr, int nc, int noffset = 0);

	// assignment operator
	void operator = (const CSRMatrix& A);

	// return row count
	int rows() const { return m_nr; }

	// return columns count
	int cols() const { return m_nc; }

	// return number of nonzeroes
	int nonzeroes() const { return (int) m_values.size(); }

	// set the value, inserting it if necessary
	void set(int i, int j, double val);

	// get a value
	double operator () (int i, int j) const;

	// see if a matrix entry was allocated
	bool isAlloc(int i, int j) const;

public:
	// matrix-vector multiplication: A.x = r
	void multv(const std::vector<double>& x, std::vector<double>& r);
	void multv(const double* x, double* r);

public:
	std::vector<double>& values() { return m_values; }
	std::vector<int>& indices() { return m_columns; }
	std::vector<int>& pointers() { return m_rowIndex; }

private:
	int		m_nr;		// number of rows
	int		m_nc;		// number of columns
	int		m_offset;	// offset (0 or 1)
	std::vector<int>	m_rowIndex;		// start of row in columns array
	std::vector<int>	m_columns;		// columns of non-zero entries
	std::vector<double>	m_values;		// values of matrix
};
