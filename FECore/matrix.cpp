// matrix.cpp: implementation of the matrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "matrix.h"
#include <assert.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void lubksb(double**a, int n, int *indx, double b[]);
void ludcmp(double**a, int n, int* indx);

vector<double> operator / (vector<double>& b, matrix& m)
{
	int n = b.size();

	vector<double> x(b);
	vector<int> indx(n);
	matrix a(m);

	ludcmp(a, n, indx);
	lubksb(a, n, indx, x);

	return x;
}

vector<double> operator * (matrix& m, vector<double>& b)
{
	int i, j;
	int NR = m.rows();
	int NC = m.columns();
	assert(NC == b.size());
	vector<double> r(NR);
	for (i=0; i<NR; ++i)
	{
		r[i] = 0.0;
		for (j=0; j<NC; ++j) r[i] += m[i][j]*b[j];
	}

	return r;
}

matrix matrix::operator * (const matrix& m)
{
	assert(m_nc == m.m_nr);
	matrix a(m_nr, m.m_nc);

	for (int i=0; i<m_nr; ++i)
	{
		for (int j=0; j<m.m_nc; ++j)
		{
			a(i,j) = 0;
			for (int k=0; k<m_nc; ++k) a(i,j) += m_pr[i][k]*m(k,j);
		}
	}

	return a;
}

matrix matrix::inverse()
{
	// make sure this is a square matrix
	assert(m_nr == m_nc);

	// make a copy of this matrix
	// since we don't want to change it
	matrix a(*this);

	// do a LU decomposition
	int n = m_nr;
	vector<int> indx(n);
	ludcmp(a, n, indx);

	// allocate the inverse matrix
	matrix ai(n, n);

	// do a backsubstituation on the columns of a
	vector<double> b(n); b.zero();
	for (int j=0; j<n; ++j)
	{
		b[j] = 1;
		lubksb(a, n, indx, b);

		for (int i=0; i<n; ++i)
		{
			ai[i][j] = b[i];
			b[i] = 0;
		}
	}

	return ai;
}
