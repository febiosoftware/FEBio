// SparseMatrix.cpp: implementation of the SparseMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SparseMatrix.h"

//-----------------------------------------------------------------------------
SparseMatrix::SparseMatrix()
{
	m_ndim = 0;
	m_nsize = 0;
	m_pd = 0;
}

//-----------------------------------------------------------------------------
void SparseMatrix::Clear()
{
	if (m_pd) delete [] m_pd; m_pd = 0;
}

//-----------------------------------------------------------------------------
void SparseMatrix::zero()
{
	memset(m_pd, 0, m_nsize*sizeof(double)); 
}

//-----------------------------------------------------------------------------
//! Print a block from a sparse matrix to file
void print(SparseMatrix& m, FILE* fp, int i0, int j0, int i1, int j1)
{
	int ndim = m.Size();
	if ((i1 < 0) || (i1 >= ndim)) i1 = ndim-1;
	if ((j1 < 0) || (j1 >= ndim)) j1 = ndim-1;

	for (int i=i0; i<=i1; ++i)
	{
		for (int j=j0; j<=j1; ++j)
		{
			fprintf(fp, "%10.3g", m.get(i,j));
		}
		fprintf(fp, "\n");
	}
}
