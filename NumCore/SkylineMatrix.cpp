#include "stdafx.h"
#include "SkylineMatrix.h"

//=============================================================================
// SkylineMatrix
//=============================================================================

//-----------------------------------------------------------------------------
SkylineMatrix::SkylineMatrix()
{
	m_ppointers = 0;
	m_pd = 0;
}

//-----------------------------------------------------------------------------
SkylineMatrix::~SkylineMatrix()
{
	Clear();
}

//-----------------------------------------------------------------------------
void SkylineMatrix::Zero()
{
	memset(m_pd, 0, m_nsize*sizeof(double)); 
}

//-----------------------------------------------------------------------------
void SkylineMatrix::Clear()
{
	if (m_pd       ) delete [] m_pd       ; m_pd = 0;
	if (m_ppointers) delete [] m_ppointers; m_ppointers = 0;

	SparseMatrix::Clear();
}

//-----------------------------------------------------------------------------
//! \todo Can I get rid of this function?
void SkylineMatrix::Create(double* pv, int* pp, int N)
{
	delete [] m_pd  ; m_pd = pv;
	delete [] m_ppointers; m_ppointers = pp;

	m_nrow = m_ncol = N;
	m_nsize = pp[N];
}

//-----------------------------------------------------------------------------
void SkylineMatrix::Create(SparseMatrixProfile& mp)
{
	int neq = mp.Rows();
	int* pointers = new int[neq + 1];

	pointers[0] = 0;
	for (int i=1; i<=neq; ++i)
	{
		vector<int>& a = mp.Column(i-1);
		int n = i - a[0];
		pointers[i] = pointers[i-1] + n;
	}

	// allocate stiffness matrix
	double* values = new double[pointers[neq]];
	if (values==0)
	{
		double falloc = (double) sizeof(double) * (double) (pointers[neq]);
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

bool SkylineMatrix::check(int i, int j)
{
	assert(false);
	return true;
}

void SkylineMatrix::add(int i, int j, double v)
{
	// only add to the upper triangular part
	if (j>=i) m_pd[ m_ppointers[j] + j-i] += v;
}

void SkylineMatrix::set(int i, int j, double v)
{
	// only add to the upper triangular part
	if (j>=i) m_pd[ m_ppointers[j] + j-i] = v;
}

double SkylineMatrix::get(int i, int j)
{
	if (i>j) { i^= j; j ^= i; i ^= j; }

	int l = m_ppointers[j+1] - m_ppointers[j];
	if (j-i < l) return m_pd[ m_ppointers[j] + j-i];
	return 0;
}

double SkylineMatrix::diag(int i)
{
	return m_pd[ m_ppointers[i] ];
}
