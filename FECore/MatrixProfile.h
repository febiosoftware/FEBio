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
//! Memory exception class
//! \todo move this to a seperate file
class MemException
{
public:
	MemException(double falloc = -1){ m_falloc = falloc; }

public:
	double m_falloc;
};

//-----------------------------------------------------------------------------
//! This class stores the profile of a sparse matrix. A profile is defined by the
//! column and row indices of the non-zero elements of a matrix. 
//! These elements are stored in a condensed format.
//! This means that for each column, an array of pairs is stored where each pair
//! identifies the start and end row index of the nonzero elements in that column.
//! The matrix profile is used to build the sparse matrix structure 
//! in an efficient way.
//! It is currently assumed that the sparse matrix is square and symmetric

class SparseMatrixProfile
{
public:
	//! Constructor. Takes the nr of equations as the input argument
	SparseMatrixProfile(int n = 0);

	//! destructor
	virtual ~SparseMatrixProfile(){}

	//! clears the matrix profile
	void clear() { m_prof.clear(); }

	//! copy constructor
	SparseMatrixProfile(SparseMatrixProfile& mp);

	//! assignment operator
	SparseMatrixProfile& operator = (SparseMatrixProfile& mp);

	//! updates the profile for an array of elements
	void UpdateProfile(vector< vector<int> >& LM, int N);

	//! returns the size of the profile. That is the nr of equations
	int size() { return m_prof.size(); }

	//! returns the non-zero row indices (in condensed format) for a column
	vector<int>& column(int i) { return m_prof[i]; }

protected:
	vector< vector<int> >	m_prof;	//!< the actual profile in condensed format
};

#endif // !defined(AFX_MATRIXPROFILE_H__F83C6F4F_AB5B_445F_AD8C_9C0CBAD26D09__INCLUDED_)
