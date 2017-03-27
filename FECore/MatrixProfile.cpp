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
SparseMatrixProfile::SparseMatrixProfile(int n)
{
	m_updateMethod = 0;

	if (n>0)
	{
		// allocate the profile array
		m_prof.resize(n);

		// initialize the profile
		for (int i=0; i<n; ++i)
		{
			vector<int>& a = m_prof[i];
			a.reserve(10);
			a.resize(2);
			a[0] = i;
			a[1] = i;
		}
	}
}

//-----------------------------------------------------------------------------
//! Copy constructor. Simply copies the profile

SparseMatrixProfile::SparseMatrixProfile(const SparseMatrixProfile& mp)
{
	m_updateMethod = mp.m_updateMethod;
	m_prof = mp.m_prof;
}

//-----------------------------------------------------------------------------
//! Assignment operator. Copies the profile.
 
SparseMatrixProfile& SparseMatrixProfile::operator =(const SparseMatrixProfile& mp)
{
	m_updateMethod = mp.m_updateMethod;
	if (m_prof.size() != mp.m_prof.size()) m_prof = mp.m_prof;
	else
	{
		for (size_t i=0; i<m_prof.size(); ++i) m_prof[i] = mp.m_prof[i];
	}

	return (*this);
}

void SparseMatrixProfile::SetUpdateMethod(SparseMatrixProfile::UpdateMethod m)
{
	m_updateMethod = (int) m;
}

void SparseMatrixProfile::UpdateProfile(vector< vector<int> >& LM, int M)
{
	switch (m_updateMethod)
	{
	case 0: UpdateProfile1(LM, M); break;
	case 1: UpdateProfile2(LM, M); break;
	default:
		assert(false);
		break;
	}
}

//-----------------------------------------------------------------------------
//! Updates the profile. The LM array contains a list of elements that contribute
//! to the sparse matrix. Each "element" defines a set of degrees of freedom that
//! are somehow connected. Each pair of dofs that are connected contributes to
//! the global stiffness matrix and therefor also to the matrix profile.

void SparseMatrixProfile::UpdateProfile1(vector< vector<int> >& LM, int M)
{
	// get the nr of equations
	int neq = (int) m_prof.size();

	// make sure there is work to do
	if (neq==0) return;

	// Count the number of elements that contribute to a certain column
	// The pval array stores this number (which I also call the valence
	// of the column)
	int* pval = new int[neq];
	if (pval == 0) throw MemException((double)(sizeof(int)*neq));

	// initial column valences to zero
	for (int i=0; i<neq; ++i) pval[i] = 0;

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
	int* pelc = new int[Ntot];
	if (pelc == 0) throw MemException(sizeof(int)*Ntot);

	int** ppelc = new int*[neq];
	if (ppelc == 0) throw MemException(sizeof(int*)*neq);

	// set the column pointers
	ppelc[0] = pelc;
	for (int i=1; i<neq; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

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
	ppelc[0] = pelc;
	for (int i=1; i<neq; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	// The pcol array is used to store the non-zero row indices
	// for a particular column.
	int* pcol = new int[neq];
	if (pcol == 0) throw MemException(sizeof(int)*neq);

	// zero the row indices for the column
	for (int j=0; j<neq; ++j) pcol[j] = -1;

	// loop over all columns
	int nold;
	for (int i=0; i<neq; ++i)
	{
		if (pval[i] > 0)
		{
			// expand column i. That is flag the non-zero rows
			// that are currently in use for this column.
			nold = 0;
			vector<int>& a = m_prof[i];
			int lmin = neq;
			for (int j=0; j<(int) a.size(); j += 2)
			{
				nold += a[j+1] - a[j] + 1;
				if (a[j] < lmin) lmin = a[j];
				for (int k=a[j]; k<=a[j+1]; ++k) pcol[k] = i;
			}
		
			// loop over all elements in the plec, flagging
			// all rows that are being used. Note that we 
			// assume a symmetric matrix so we only store the
			// profile of the upper triangular matrix.
			int n = 0;
			for (int j=0; j<pval[i]; ++j)
			{
				int iel = (ppelc[i])[j];
				int* lm = &(LM[iel])[0];
				int N = LM[iel].size();
				for (int k=0; k<N; ++k)
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

				int l = lmin;
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

//-----------------------------------------------------------------------------
// extract the matrix profile of a block
SparseMatrixProfile SparseMatrixProfile::GetBlockProfile(int nrow0, int ncol0, int nrow1, int ncol1) const
{
	// This will store the block profile
	SparseMatrixProfile bMP;

	// number of columns in block
	int NC = ncol1 - ncol0 + 1;
	assert(NC > 0);
	bMP.m_prof.resize(NC);

	for (int j=0; j<NC; ++j)
	{
		const vector<int>& sj = m_prof[ncol0+j];
		vector<int>& dj = bMP.m_prof[j];
		int nr = sj.size();
		for (int i=0; i<nr; i+=2)
		{
			int n0 = sj[2*i  ];
			int n1 = sj[2*i+1];
			assert(n0<=n1);

			if ((n1 >= nrow0)&&(n0 <= nrow1))
			{
				if (n0 < nrow0) n0 = nrow0;
				if (n1 > nrow1) n1 = nrow1;

				dj.push_back(n0);
				dj.push_back(n1);
			}
		}
	}

	return bMP;
}

//-----------------------------------------------------------------------------
// this sort function is defined in qsort.cpp
void qsort(int n, int* arr, int* indx);

//-----------------------------------------------------------------------------
void SparseMatrixProfile::UpdateProfile2(vector< vector<int> >& LM, int lmsize)
{
	// get the nr of equations
	int neq = m_prof.size();

	// make sure there is work to do
	if (neq == 0) return;

	// loop over all elements
	for (int i = 0; i<lmsize; ++i)
	{
		// get the element indices
		vector<int>& lm = LM[i];

		// get the size of this element's array
		const int M = (int) lm.size();
		if (M == 0) continue;

		// sort the indices
		vector<int> P(M);
		qsort(M, &lm[0], &P[0]);

		// loop over all columns
		for (int j=0; j<M; ++j)
		{
			int ncol = lm[P[j]];
			if (ncol >= 0)
			{
				// get the column data and its size
				vector<int>& a = m_prof[ncol];
				int N = (int)a.size();

				// loop over all indices
				int n = 0;
				for (int m = 0; m<M; ++m)
				{
					// get the next row index
					int nrow = lm[P[m]];

					// we only store upper triangular profile
					if ((nrow >= 0) && (nrow < ncol))
					{
						bool bdone = false;
						do
						{
							int n0 = a[n];	// start row
							int n1 = a[n + 1];	// end row
							int n2 = (n<N - 2 ? a[n + 2] : -1);	// start row of next pair
							assert(n1 >= n0);
							assert((n2 == -1) || (n2 > n1 + 1));

							// see if the row is already inserted
							if ((nrow < n0) || (nrow > n1))
							{
								// check for merge
								if ((nrow == n1 + 1) && (nrow == n2 - 1))
								{
									// we must merge the two pairs
									assert(n<N - 2);
									a[n + 1] = a[n + 3];
									a.erase(a.begin() + n + 2, a.begin() + n + 4);
									N -= 2;
									bdone = true;
								}
								else if (nrow == n0 - 1)
								{
									// see if we can simply grow the interval
									a[n] = nrow;
									bdone = true;
								}
								else if (nrow == n1 + 1)
								{
									a[n + 1] = nrow;
									bdone = true;
								}
								else if (nrow < n0)
								{
									// see if we need to insert the row
									N += 2;
									a.insert(a.begin() + n, 2, nrow);
									bdone = true;
								}
								else n += 2;
							}
							else bdone = true;
						} while ((bdone == false) && (n<N));

						if (bdone == false)
						{
							a.push_back(nrow);
							a.push_back(nrow);
							N += 2;
						}
					}
				}
			}
		}
	}
}
