#include "stdafx.h"
#include "SkylineMatrix.h"

//=============================================================================
// SkylineMatrix
//=============================================================================

//-----------------------------------------------------------------------------
SkylineMatrix::SkylineMatrix()
{
	m_ppointers = 0;
}

//-----------------------------------------------------------------------------
SkylineMatrix::~SkylineMatrix()
{
	delete [] m_pd;
	delete [] m_ppointers;
}

//-----------------------------------------------------------------------------
//! \todo Can I get rid of this function?
void SkylineMatrix::Create(double* pv, int* pp, int N)
{
	delete [] m_pd  ; m_pd = pv;
	delete [] m_ppointers; m_ppointers = pp;

	m_ndim  = N;
	m_nsize = pp[N];
}

//-----------------------------------------------------------------------------
void SkylineMatrix::Create(SparseMatrixProfile& mp)
{
	int i, n;

	int neq = mp.size();
	int* pointers = new int[neq + 1];
	if (pointers == 0) throw MemException(sizeof(int)*(neq+1));

	pointers[0] = 0;
	for (i=1; i<=neq; ++i)
	{
		vector<int>& a = mp.column(i-1);
		n = i - a[0];
		pointers[i] = pointers[i-1] + n;
	}

	// allocate stiffness matrix
	double* values = new double[pointers[neq]];
	if (values==0)
	{
		double falloc = (double) sizeof(double) * (double) (pointers[neq]);
		throw MemException(falloc);
	}

	// create the matrix
	Create(values, pointers, neq);
}


//-----------------------------------------------------------------------------
//! This function assembles the local stiffness matrix
//! into the global stiffness matrix which is in skyline format
//!
void SkylineMatrix::Assemble(matrix& ke, vector<int>& LM)
{
	int i, j, I, J;

	const int N = ke.rows();

	double* pv = values();
	int* pi = pointers();

	for (i=0; i<N; ++i)
	{
		I = LM[i];

		if (I>=0)
		{
			for (j=0; j<N; ++j)
			{
				J = LM[j];

				// only add values to upper-diagonal part of stiffness matrix
				if (J>=I)
				{
					pv[ pi[J] + J - I] += ke[i][j];
				}
			}
		}
	}
}


//-----------------------------------------------------------------------------
void SkylineMatrix::Assemble(matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	int i, j, I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	double* pv = values();
	int* pi = pointers();

	for (i=0; i<N; ++i)
	{
		I = LMi[i];

		if (I>=0)
		{
			for (j=0; j<M; ++j)
			{
				J = LMj[j];

				// only add values to upper-diagonal part of stiffness matrix
				if (J>=I)
				{
					pv[ pi[J] + J - I] += ke[i][j];
				}
			}
		}
	}
}
