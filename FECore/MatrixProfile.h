// MatrixProfile.h: interface for the MatrixProfile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIXPROFILE_H__F83C6F4F_AB5B_445F_AD8C_9C0CBAD26D09__INCLUDED_)
#define AFX_MATRIXPROFILE_H__F83C6F4F_AB5B_445F_AD8C_9C0CBAD26D09__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

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
		vector<RowEntry>	m_data;	// the column profile data
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
	void UpdateProfile(vector< vector<int> >& LM, int N);

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
	vector<ColumnProfile>	m_prof;	//!< the actual profile in condensed format
};

#endif // !defined(AFX_MATRIXPROFILE_H__F83C6F4F_AB5B_445F_AD8C_9C0CBAD26D09__INCLUDED_)
