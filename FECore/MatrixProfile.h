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
#include "vector.h"

//-----------------------------------------------------------------------------
//! This class stores the profile of a sparse matrix. A profile is defined by the
//! column and row indices of the non-zero elements of a matrix. 
//! These elements are stored in a condensed format.
//! This means that for each column, an array of pairs is stored where each pair
//! identifies the start and end row index of the nonzero elements in that column.
//! The matrix profile is used to build the sparse matrix structure 
//! in an efficient way.

class FECORE_API SparseMatrixProfile
{
public:
	struct RowEntry
	{
		int start, end;
	};

	class ColumnProfile
	{
	public:
		ColumnProfile() {}

		ColumnProfile(const ColumnProfile& a);

		// get the number of row entries
		int size() const { return (int) m_data.size(); }

		// access
		RowEntry& operator [] (int i) { return m_data[i]; }
		const RowEntry& operator [] (int i) const { return m_data[i]; }

		// make room
		void clear() { m_data.clear(); }

		// reserve some storage
		void reserve(int n)
		{
			m_data.reserve(n);
		}

		// add to the end
		void push_back(int n0, int n1)
		{
			RowEntry re = { n0, n1 };
			m_data.push_back(re);
		}

		void push_front(int n0, int n1)
		{
			RowEntry re = {n0, n1};
			m_data.insert(m_data.begin(), re);
		}

		// add row index to column profile
		void insertRow(int row);

	private:
		std::vector<RowEntry>	m_data;	// the column profile data
	};

public:
	//! Constructor. Takes the nr of equations as the input argument
	SparseMatrixProfile(int nrow = 0, int ncol = 0);

	//! allocate storage for profile
	void Create(int nrow, int ncol);

	//! copy constructor
	SparseMatrixProfile(const SparseMatrixProfile& mp);

	//! assignment operator
	SparseMatrixProfile& operator = (const SparseMatrixProfile& mp);

	//! Create the profile of a diagonal matrix
	void CreateDiagonal();

	//! clears the matrix profile
	void Clear();

	//! updates the profile for an array of elements
	void UpdateProfile(std::vector< std::vector<int> >& LM, int N);

	//! inserts an entry into the profile (This is an expensive operation!)
	void Insert(int i, int j);

	//! returns the number of rows
	int Rows() const { return m_nrow; }

	//! returns the number of columns
	int Columns() const { return m_ncol; }

	//! returns the non-zero row indices (in condensed format) for a column
	ColumnProfile& Column(int i) { return m_prof[i]; }

	// Extracts a block profile
	SparseMatrixProfile GetBlockProfile(int nrow0, int ncol0, int nrow1, int ncol1) const;

private:
	int	m_nrow, m_ncol;				//!< dimensions of matrix
	std::vector<ColumnProfile>	m_prof;	//!< the actual profile in condensed format
};
