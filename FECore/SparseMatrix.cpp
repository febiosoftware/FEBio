// SparseMatrix.cpp: implementation of the SparseMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SparseMatrix.h"

//////////////////////////////////////////////////////////////////////
// SparseMatrix
//////////////////////////////////////////////////////////////////////

SparseMatrix::SparseMatrix()
{
	m_ndim = 0;
	m_nsize = 0;
	m_pd = 0;
}

SparseMatrix::~SparseMatrix()
{

}

void SparseMatrix::print(FILE* fp, int i0, int j0, int i1, int j1)
{
	if ((i1 < 0) || (i1 >= m_ndim)) i1 = m_ndim-1;
	if ((j1 < 0) || (j1 >= m_ndim)) j1 = m_ndim-1;

	for (int i=i0; i<=i1; ++i)
	{
		for (int j=j0; j<=j1; ++j)
		{
			fprintf(fp, "%10.3g", get(i,j));
		}
		fprintf(fp, "\n");
	}
}

//////////////////////////////////////////////////////////////////////
// FullMatrix
//////////////////////////////////////////////////////////////////////

FullMatrix::FullMatrix()
{
	m_pr = 0;
}

FullMatrix::~FullMatrix()
{
	delete [] m_pd; m_pd = 0;
	delete [] m_pr; m_pr = 0;
}

void FullMatrix::Create(int N)
{
	if (N != m_ndim)
	{
		if (m_pd) delete [] m_pd;
		if (m_pr) delete [] m_pr;

		m_pd = new double[N*N];
		m_pr = new double*[N];

		for (int i=0; i<N; ++i) m_pr[i] = m_pd + i*N;

		m_ndim = N;
		m_nsize = N*N;
	}
}

//////////////////////////////////////////////////////////////////////
// SkylineMatrix
//////////////////////////////////////////////////////////////////////

SkylineMatrix::SkylineMatrix()
{
	m_ppointers = 0;
}

SkylineMatrix::~SkylineMatrix()
{
	delete [] m_pd;
	delete [] m_ppointers;
}

void SkylineMatrix::Create(double* pv, int* pp, int N)
{
	delete [] m_pd  ; m_pd = pv;
	delete [] m_ppointers; m_ppointers = pp;

	m_ndim  = N;
	m_nsize = pp[N];
}

//////////////////////////////////////////////////////////////////////
// CompactMatrix
//////////////////////////////////////////////////////////////////////

CompactMatrix::CompactMatrix(int offset)
{
	m_pindices = 0;
	m_ppointers = 0;
	m_offset = offset;
}

void CompactMatrix::Create(int N, int nz, double* pv, int* pi, int* pp)
{
	if (m_pd  ) delete [] m_pd; m_pd = pv;
	if (m_pindices ) delete [] m_pindices; m_pindices = pi;
	if (m_ppointers) delete [] m_ppointers; m_ppointers = pp;

	m_ndim  = N;
	m_nsize = nz;
}


bool CompactMatrix::print_hb()
{
	int isize, dsize;

	isize = sizeof(int);
	dsize = sizeof(double);

	FILE* fout = fopen("hb_matrix.out", "wb");
	if (fout == 0) { fprintf(stderr, "Failed creating output file."); return false; }

	fwrite(&m_ndim, isize, 1, fout);
	fwrite(&m_nsize, isize, 1, fout);
	fwrite(m_ppointers, isize, m_ndim+1, fout);
	fwrite(m_pindices, isize, m_nsize, fout);
	fwrite(m_pd, dsize, m_nsize, fout);

	fclose(fout);
	return true;

}

//////////////////////////////////////////////////////////////////////
// CompactSymmMatrix
//////////////////////////////////////////////////////////////////////

CompactSymmMatrix::CompactSymmMatrix(int offset) : CompactMatrix(offset) {}

void CompactSymmMatrix::mult_vector(const vector<double>& x, vector<double>& r)
{
	int j, i, n;
	int N = x.size();
	assert(N == Size());

	double* pv, rj;
	int* pi;

	r.zero();
	// loop over all columns
	for (j=0; j<N; ++j)
	{
		pv = m_pd  + m_ppointers[j];
		pi = m_pindices + m_ppointers[j];
		n = m_ppointers[j+1] - m_ppointers[j];

		// add off-diagonal elements
		for (i=1; i<n-7; i+=8)
		{
			// add lower triangular element
			r[pi[i  ]] += pv[i  ]*x[j];
			r[pi[i+1]] += pv[i+1]*x[j];
			r[pi[i+2]] += pv[i+2]*x[j];
			r[pi[i+3]] += pv[i+3]*x[j];
			r[pi[i+4]] += pv[i+4]*x[j];
			r[pi[i+5]] += pv[i+5]*x[j];
			r[pi[i+6]] += pv[i+6]*x[j];
			r[pi[i+7]] += pv[i+7]*x[j];
		}
		for (i=0; i<(n-1)%8; ++i)
			r[pi[n-1-i]] += pv[n-1-i]*x[j];

		rj = pv[0]*x[j]; // add diagonal element
		for (i=1; i<n-7; i+=8)
		{
			// add upper triangular element
			rj += pv[i  ]*x[pi[i  ]];
			rj += pv[i+1]*x[pi[i+1]];
			rj += pv[i+2]*x[pi[i+2]];
			rj += pv[i+3]*x[pi[i+3]];
			rj += pv[i+4]*x[pi[i+4]];
			rj += pv[i+5]*x[pi[i+5]];
			rj += pv[i+6]*x[pi[i+6]];
			rj += pv[i+7]*x[pi[i+7]];
		}
		for (i=0; i<(n-1)%8; ++i)
			rj += pv[n-1-i]*x[pi[n-1-i]];
		r[j] += rj;
	}
}

CompactUnSymmMatrix::CompactUnSymmMatrix(int offset, bool row_based) : CompactMatrix(offset)
	{
		m_brow_based = row_based;
	}
