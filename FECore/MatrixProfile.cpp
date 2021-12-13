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



#include "stdafx.h"
#include "MatrixProfile.h"
#include <assert.h>
using namespace std;

SparseMatrixProfile::ColumnProfile::ColumnProfile(const SparseMatrixProfile::ColumnProfile& a)
{
	m_data = a.m_data;
}

void SparseMatrixProfile::ColumnProfile::insertRow(int row)
{
	// first, check if empty
	if (m_data.empty()) { push_back(row, row); return; }

	int N = size();

	// check some easy cases first
	if (row + 1 <  m_data[0].start) { push_front(row, row); return; }
	if (row + 1 == m_data[0].start) { m_data[0].start--; return; }

	if (row - 1 >  m_data[N-1].end) { push_back(row, row); return; }
	if (row - 1 == m_data[N-1].end) { m_data[N-1].end++; return; }

	// general case, find via bisection
	int N0 = 0, N1 = N-1;
	int n = N/2, m = 0;
	while (true)
	{
		RowEntry& rn = m_data[n];

		// see if row falls inside the interval
		if ((row >= rn.start) && (row <= rn.end))
		{
			// no need to do anything
			return;
		}

		if (row < rn.start)
		{
			assert(n>0); // this should always be the case due to the easy case handling above

			// get the previous entry
			RowEntry& r0 = m_data[n-1];

			if (row > r0.end)
			{
				if (row + 1 == rn.start)
				{
					if (r0.end == row - 1)
					{
						// merge entries
						r0.end = rn.end;
						m_data.erase(m_data.begin() + n);
						return;
					}
					else 
					{
						rn.start--;
						return;
					}
				}
				else if (row - 1 == r0.end)
				{
					r0.end++;
					return;
				}
				else
				{
					RowEntry re = {row, row};
					m_data.insert(m_data.begin() + n, re);
					return;
				}
			}
			else
			{
				N1 = n;
				n = (N0 + N1) / 2;
			}
		}
		else
		{
			assert(row > rn.end);

			// get the next entry
			RowEntry& r1 = m_data[n + 1];

			if (row < r1.start)
			{
				if (row - 1 == rn.end)
				{
					if (r1.start == row + 1)
					{
						// merge entries
						r1.start = rn.start;
						m_data.erase(m_data.begin() + n);
						return;
					}
					else
					{
						rn.end++;
						return;
					}
				}
				else if (row + 1 == r1.start)
				{
					r1.start--;
					return;
				}
				else
				{
					RowEntry re = { row, row };
					m_data.insert(m_data.begin() + n + 1, re);
					return;
				}
			}
			else
			{
				N0 = n;
				n = (N0 + N1 + 1) / 2;
			}
		}
		++m;
		assert(m <= N);
	}
}

//-----------------------------------------------------------------------------
//! MatrixProfile constructor. Takes the nr of equations as input argument.
//! If n is larger than zero a default profile is constructor for a diagonal
//! matrix.
SparseMatrixProfile::SparseMatrixProfile(int nrow, int ncol)
{
	m_nrow = nrow;
	m_ncol = ncol;

	// allocate storage profile
	if (ncol > 0) 
	{
		int nres = (m_ncol < 100 ? m_ncol : 100);
		m_prof.resize(ncol);
		for (int i=0; i<ncol; ++i) m_prof[i].reserve(nres);
	}
}

//-----------------------------------------------------------------------------
//! allocate storage for profile
void SparseMatrixProfile::Create(int nrow, int ncol)
{
	m_nrow = nrow;
	m_ncol = ncol;

	int nres = (m_ncol < 100 ? m_ncol : 100);
	m_prof.resize(ncol);
	for (int i = 0; i<ncol; ++i) m_prof[i].reserve(nres);
}

//-----------------------------------------------------------------------------
//! Copy constructor. Simply copies the profile

SparseMatrixProfile::SparseMatrixProfile(const SparseMatrixProfile& mp)
{
	m_nrow = mp.m_nrow;
	m_ncol = mp.m_ncol;
	m_prof = mp.m_prof;
}

//-----------------------------------------------------------------------------
//! Assignment operator. Copies the profile.
 
SparseMatrixProfile& SparseMatrixProfile::operator =(const SparseMatrixProfile& mp)
{
	m_nrow = mp.m_nrow;
	m_ncol = mp.m_ncol;
	m_prof = mp.m_prof;

	return (*this);
}

//-----------------------------------------------------------------------------
//! Create the profile of a diagonal matrix
void SparseMatrixProfile::CreateDiagonal()
{
	int n = min(m_nrow, m_ncol);

	// initialize the profile to a diagonal matrix
	for (int i = 0; i<n; ++i)
	{
		ColumnProfile& a = m_prof[i];
		a.insertRow(i);
	}
}

//-----------------------------------------------------------------------------
void SparseMatrixProfile::Clear()
{ 
	m_prof.clear(); 
}

//-----------------------------------------------------------------------------
//! Updates the profile. The LM array contains a list of elements that contribute
//! to the sparse matrix. Each "element" defines a set of degrees of freedom that
//! are somehow connected. Each pair of dofs that are connected contributes to
//! the global stiffness matrix and therefor also to the matrix profile.
void SparseMatrixProfile::UpdateProfile(vector< vector<int> >& LM, int M)
{
	// get the dimensions of the matrix
	int nr = m_nrow;
	int nc = m_ncol;

	// make sure there is work to do
	if (nr*nc == 0) return;

	// Count the number of elements that contribute to a certain column
	// The pval array stores this number (which I also call the valence
	// of the column)
	vector<int> pval(nc, 0);

	// fill the valence array
	int Ntot = 0;
	for (int i = 0; i<M; ++i)
	{
		int* lm = &(LM[i])[0];
		int N = (int)LM[i].size();
		Ntot += N;
		for (int j = 0; j<N; ++j)
		{
			if (lm[j] >= 0) pval[lm[j]]++;
		}
	}

	// create a "compact" 2D array that stores for each column the element
	// numbers that contribute to that column. The compact array consists
	// of two arrays. The first one (pelc) contains all element numbers, sorted
	// by column. The second array stores for each column a pointer to the first
	// element in the pelc array that contributes to that column.
	vector<int> pelc(Ntot);
	vector<int*> ppelc(nc);

	// set the column pointers
	ppelc[0] = &pelc[0];
	for (int i = 1; i<nc; ++i) ppelc[i] = ppelc[i - 1] + pval[i - 1];

	// fill the pelc array
	for (int i = 0; i<M; ++i)
	{
		int* lm = &(LM[i])[0];
		int N = (int)LM[i].size();
		for (int j = 0; j<N; ++j)
		{
			if (lm[j] >= 0) *(ppelc[lm[j]])++ = i;
		}
	}

	// reset pelc pointers
	ppelc[0] = &pelc[0];
	for (int i = 1; i<nc; ++i) ppelc[i] = ppelc[i - 1] + pval[i - 1];

	// loop over all columns
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<nc; ++i)
	{
		if (pval[i] > 0)
		{
			// get the column
			ColumnProfile& a = m_prof[i];

			// loop over all elements in the plec
			for (int j = 0; j<pval[i]; ++j)
			{
				int iel = (ppelc[i])[j];
				int* lm = &(LM[iel])[0];
				int N = (int)LM[iel].size();
				for (int k = 0; k<N; ++k)
				{
					if (lm[k] >= 0)
					{ 
						a.insertRow(lm[k]);
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! inserts an entry into the profile
void SparseMatrixProfile::Insert(int i, int j)
{
	ColumnProfile& a = m_prof[j];
	a.insertRow(i);
}

//-----------------------------------------------------------------------------
// extract the matrix profile of a block
SparseMatrixProfile SparseMatrixProfile::GetBlockProfile(int nrow0, int ncol0, int nrow1, int ncol1) const
{
	int nrows = nrow1 - nrow0 + 1;
	int ncols = ncol1 - ncol0 + 1;
	assert(nrows > 0);
	assert(ncols > 0);

	// This will store the block profile
	SparseMatrixProfile bMP(nrows, ncols);

	for (int j=0; j<ncols; ++j)
	{
		const ColumnProfile& sj = m_prof[ncol0+j];
		ColumnProfile& dj = bMP.m_prof[j];
		int nr = sj.size();
		for (int i=0; i<nr; i++)
		{
			const RowEntry& ri = sj[i];

			int n0 = ri.start;
			int n1 = ri.end;

			if ((n1 >= nrow0)&&(n0 <= nrow1))
			{
				if (n0 < nrow0) n0 = nrow0;
				if (n1 > nrow1) n1 = nrow1;

				n0 -= nrow0;
				n1 -= nrow0;

				dj.push_back(n0, n1);
			}
		}
	}

	return bMP;
}
