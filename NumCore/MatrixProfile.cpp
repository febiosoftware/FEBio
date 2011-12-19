// MatrixProfile.cpp: implementation of the MatrixProfile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MatrixProfile.h"

//-----------------------------------------------------------------------------
//! MatrixProfile constructor. Takes the nr of equations as input argument.
//! If n is larger than zero a default profile is constructor for a diagonal
//! matrix.
NumCore::SparseMatrixProfile::SparseMatrixProfile(int n)
{
	if (n>0)
	{
		// allocate the profile array
		m_prof.resize(n);

		// initialize the profile
		for (int i=0; i<n; ++i)
		{
			vector<int>& a = m_prof[i];
			a.resize(2);
			a[0] = i;
			a[1] = i;
		}
	}
}

//-----------------------------------------------------------------------------
//! Copy constructor. Simply copies the profile

NumCore::SparseMatrixProfile::SparseMatrixProfile(SparseMatrixProfile& mp)
{
	m_prof = mp.m_prof;
}

//-----------------------------------------------------------------------------
//! Assignment operator. Copies the profile.
 
NumCore::SparseMatrixProfile& NumCore::SparseMatrixProfile::operator =(SparseMatrixProfile& mp)
{
	if (m_prof.size() != mp.m_prof.size()) m_prof = mp.m_prof;
	else
	{
		for (size_t i=0; i<m_prof.size(); ++i) m_prof[i] = mp.m_prof[i];
	}

	return (*this);
}

//-----------------------------------------------------------------------------
//! Updates the profile. The LM array contains a list of elements that contribute
//! to the sparse matrix. Each "element" defines a set of degrees of freedom that
//! are somehow connected. Each pair of dofs that are connected contributes to
//! the global stiffness matrix and therefor also to the matrix profile.

void NumCore::SparseMatrixProfile::UpdateProfile(vector< vector<int> >& LM, int M)
{
	int i, j, k, l, iel, n, *lm, N, Ntot = 0;

	// get the nr of equations
	int neq = m_prof.size();

	// Count the number of elements that contribute to a certain column
	// The pval array stores this number (which I also call the valence
	// of the column)
	int* pval = new int[neq];
	if (pval == 0) throw MemException(sizeof(int)*neq);

	// initial column valences to zero
	for (i=0; i<neq; ++i) pval[i] = 0;

	// fill the valence array
	for (i=0; i<M; ++i)
	{
		lm = &(LM[i])[0];
		N = LM[i].size();
		Ntot += N;
		for (j=0; j<N; ++j)
		{
			if (lm[j] >= 0) pval[lm[j]]++; 
		}
	}

	// create a "compact" 2D array that stores for each column the element
	// numbers that contribute to that column. The compact array consists
	// of two arrays. The first one (pelc) contains all element numbers, sorted
	// by column. The second array stores for each column a pointer to the first
	// element in the pelc array that contributes to that column.
	int* pelc = new int[Ntot];
	if (pelc == 0) throw MemException(sizeof(int)*Ntot);

	int** ppelc = new int*[neq];
	if (ppelc == 0) throw MemException(sizeof(int*)*neq);

	// set the column pointers
	ppelc[0] = pelc;
	for (i=1; i<neq; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	// fill the pelc array
	for (i=0; i<M; ++i)
	{
		lm = &(LM[i])[0];
		N = LM[i].size();
		for (j=0; j<N; ++j)
		{
			if (lm[j] >= 0) *(ppelc[lm[j]])++ = i;
		}
	}

	// reset pelc pointers
	ppelc[0] = pelc;
	for (i=1; i<neq; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	// The pcol array is used to store the non-zero row indices
	// for a particular column.
	int* pcol = new int[neq];
	if (pcol == 0) throw MemException(sizeof(int)*neq);

	// zero the row indices for the column
	for (j=0; j<neq; ++j) pcol[j] = -1;

	// loop over all columns
	int nold;
	for (i=0; i<neq; ++i)
	{
		if (pval[i] > 0)
		{

			// expand column i. That is flag the non-zero rows
			// that are currently in use for this column.
			nold = 0;
			vector<int>& a = m_prof[i];
			int lmin = neq;
			for (j=0; j<(int) a.size(); j += 2)
			{
				nold += a[j+1] - a[j] + 1;
				if (a[j] < lmin) lmin = a[j];
				for (k=a[j]; k<=a[j+1]; ++k) pcol[k] = i;
			}
		
			// loop over all elements in the plec, flagging
			// all rows that are being used. Note that we 
			// assume a symmetric matrix so we only store the
			// profile of the upper triangular matrix.
			n = 0;
			for (j=0; j<pval[i]; ++j)
			{
				iel = (ppelc[i])[j];
				lm = &(LM[iel])[0];
				N = LM[iel].size();
				for (k=0; k<N; ++k)
				{
					if ((lm[k] >= 0) && (lm[k] < i) && (pcol[ lm[k] ] != i)) { ++n; pcol[ lm[k] ] = i; if (lm[k]<lmin) lmin = lm[k]; }
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

				l = lmin;
				do
				{
					// find a start row index
					while ((pcol[l] != i) && (l<=i)) ++l;

					// find the corresponding end row index
					if (l<=i)
					{
						a.push_back(l);
						while ((l<i) && (pcol[l+1] == i)) ++l;
						a.push_back(l);
						++l;
					}
				}
				while (l<i);
			}
		}
	}

	// cleanup temp data structures
	delete [] pcol;
	delete [] ppelc;
	delete [] pelc;
	delete [] pval;
}
