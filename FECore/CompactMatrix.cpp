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
#include "CompactMatrix.h"
#include <assert.h>

//=============================================================================
// CompactMatrix
//=============================================================================

//-----------------------------------------------------------------------------
CompactMatrix::CompactMatrix(int offset)
{
	m_pd = 0;
	m_pindices = 0;
	m_ppointers = 0;
	m_offset = offset;

	m_bdel = false;
}


//-----------------------------------------------------------------------------
CompactMatrix::~CompactMatrix()
{ 
	Clear(); 
}

//-----------------------------------------------------------------------------
void CompactMatrix::Zero()
{
	memset(m_pd, 0, m_nsize*sizeof(double)); 
}

//-----------------------------------------------------------------------------
void CompactMatrix::Clear()
{
	if (m_bdel)
	{
		if (m_pd) delete [] m_pd;
		if (m_pindices) delete [] m_pindices;
		if (m_ppointers) delete [] m_ppointers;
	}

	m_pd = 0;
	m_pindices = 0;
	m_ppointers = 0;

	SparseMatrix::Clear();
}

//-----------------------------------------------------------------------------
void CompactMatrix::alloc(int nr, int nc, int nz, double* pv, int* pi, int* pp, bool bdel)
{
	Clear();
	m_pd = pv;
	m_pindices = pi;
	m_ppointers = pp;

	m_bdel = bdel;

	m_nrow = nr;
	m_ncol = nc;
	m_nsize = nz;

	int nn = (isRowBased() ? nr : nc) + 1;
}

//-----------------------------------------------------------------------------
//! calculate bandwidth of matrix
int CompactMatrix::bandWidth()
{
	int NR = Rows();
	int NC = Columns();
	int kmax = 0;
	if (isRowBased())
	{
		for (int i = 0; i < NR; ++i)
		{
			int* pi = m_pindices + (m_ppointers[i] - m_offset);
			int n = m_ppointers[i + 1] - m_ppointers[i];
			for (int k = 0; k < n; ++k)
			{
				int j = pi[k] - m_offset;
				if (abs(j - i) > kmax) kmax = abs(j - i);
			}
		}
	}
	else
	{
		for (int j = 0; j < NC; ++j)
		{
			int* pj = m_pindices + (m_ppointers[j] - m_offset);
			int n = m_ppointers[j + 1] - m_ppointers[j];
			for (int k = 0; k < n; ++k)
			{
				int i = pj[k] - m_offset;
				if (abs(j - i) > kmax) kmax = abs(j - i);
			}
		}
	}

	return kmax;
}
