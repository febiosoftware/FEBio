// MatrixProfile.cpp: implementation of the MatrixProfile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MatrixProfile.h"
#include <assert.h>

//-----------------------------------------------------------------------------
//! MatrixProfile constructor. Takes the nr of equations as input argument.
//! If n is larger than zero a default profile is constructor for a diagonal
//! matrix.
SparseMatrixProfile::SparseMatrixProfile(int nrow, int ncol)
{
	m_nrow = nrow;
	m_ncol = ncol;

	// allocate storage profile
	if (ncol > 0) m_prof.resize(ncol);
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
		vector<int>& a = m_prof[i];
		a.reserve(10);
		a.resize(2);
		a[0] = i;
		a[1] = i;
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
	if (nr*nc==0) return;

	// Count the number of elements that contribute to a certain column
	// The pval array stores this number (which I also call the valence
	// of the column)
	vector<int> pval(nc, 0);

	// fill the valence array
	int Ntot = 0;
	for (int i=0; i<M; ++i)
	{
		int* lm = &(LM[i])[0];
		int N = LM[i].size();
		Ntot += N;
		for (int j=0; j<N; ++j)
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
	for (int i=1; i<nc; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	// fill the pelc array
	for (int i=0; i<M; ++i)
	{
		int* lm = &(LM[i])[0];
		int N = LM[i].size();
		for (int j=0; j<N; ++j)
		{
			if (lm[j] >= 0) *(ppelc[lm[j]])++ = i;
		}
	}

	// reset pelc pointers
	ppelc[0] = &pelc[0];
	for (int i=1; i<nc; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	// The pcol array is used to store the non-zero row indices
	// for a particular column.
	vector<int> pcol(nr, -1);

	// loop over all columns
	#pragma omp parallel for private(pcol) schedule(dynamic)
	for (int i=0; i<nc; ++i)
	{
		if (pval[i] > 0)
		{
			// expand column i. That is flag the non-zero rows
			// that are currently in use for this column.
			int nold = 0;
			vector<int>& a = m_prof[i];
			int lmin = nr;
			for (int j=0; j<(int) a.size(); j += 2)
			{
				nold += a[j+1] - a[j] + 1;
				if (a[j] < lmin) lmin = a[j];
				for (int k=a[j]; k<=a[j+1]; ++k) pcol[k] = i;
			}
		
			// loop over all elements in the plec, flagging
			// all rows that are being used. 
			int n = 0;
			for (int j=0; j<pval[i]; ++j)
			{
				int iel = (ppelc[i])[j];
				int* lm = &(LM[iel])[0];
				int N = LM[iel].size();
				for (int k=0; k<N; ++k)
				{
					if ((lm[k] >= 0) && (pcol[ lm[k] ] != i)) { ++n; pcol[ lm[k] ] = i; if (lm[k]<lmin) lmin = lm[k]; }
				}
			}

			// repack the column. That is, pack the non-zero row indices
			// in a condensed format. The condensed format consists of an
			// array of pairs where each pair marks the start and end row
			// of a block of non-zero matrix elements.
			if (n > 0)
			{
				// initialize the packed array to zero
				a.clear();

				int l = lmin;
				do
				{
					// find a start row index
					while ((l<nr) && (pcol[l] != i)) ++l;

					// find the corresponding end row index
					if (l<nr)
					{
						a.push_back(l);
						while ((l<nr-1) && (pcol[l+1] == i)) ++l;
						a.push_back(l);
						++l;
					}
				}
				while (l<nr);
			}
		}
	}
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
		const vector<int>& sj = m_prof[ncol0+j];
		vector<int>& dj = bMP.m_prof[j];
		int nr = sj.size();
		for (int i=0; i<nr; i+=2)
		{
			int n0 = sj[i  ];
			int n1 = sj[i+1];
			assert(n0<=n1);

			if ((n1 >= nrow0)&&(n0 <= nrow1))
			{
				if (n0 < nrow0) n0 = nrow0;
				if (n1 > nrow1) n1 = nrow1;

				n0 -= nrow0;
				n1 -= nrow0;

				dj.push_back(n0);
				dj.push_back(n1);
			}
		}
	}

	return bMP;
}
