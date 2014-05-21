// SparseMatrix.cpp: implementation of the SparseMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SparseMatrix.h"
#include <memory.h>

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
