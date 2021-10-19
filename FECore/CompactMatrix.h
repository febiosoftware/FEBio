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
#include "CSRMatrix.h"

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format.

//! This is the base class for the symmetric and unsymmetric classes

class FECORE_API CompactMatrix : public SparseMatrix
{
public:
	//! constructor
	CompactMatrix( int offset );

	//! destructor
	virtual ~CompactMatrix();

	//! zero matrix elements
	void Zero() override;

	//! Clear 
	void Clear() override;

public:
	//! Pointer to matrix values
	double* Values  () override { return m_pd;   }

	//! Pointer to matrix indices
	int*    Indices() override { return m_pindices; }

	//! pointer to matrix row pointers
	int*    Pointers() override { return m_ppointers; }

	//! return the index offset (is 0 or 1)
	int     Offset() const override { return m_offset; }

public:
	//! Create the matrix
	void alloc(int nr, int nc, int nz, double* pv, int *pi, int* pp, bool bdel = true);

	//! is the matrix symmetric or not
	virtual bool isSymmetric() = 0;

	//! is this a row-based format or not
	virtual bool isRowBased() = 0;

public:
	//! Calculate the infinity norm
	virtual double infNorm() const = 0;

	//! calculate the one norm
	virtual double oneNorm() const = 0;

	//! calculate bandwidth of matrix
	int bandWidth();

protected:
	double*	m_pd;			//!< matrix values
	int*	m_pindices;		//!< indices
	int*	m_ppointers;	//!< pointers
	int		m_offset;		//!< adjust array indices for fortran arrays
	bool	m_bdel;			//!< delete data arrays in destructor

protected:
	std::vector<int>	P;
};
