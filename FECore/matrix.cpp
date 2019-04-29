/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "matrix.h"
#include <assert.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// These functions are defined in colsol.cpp
void lubksb(double**a, int n, int *indx, double b[]);
void ludcmp(double**a, int n, int* indx);

//-----------------------------------------------------------------------------
// These functions are defined in svd.cpp
void svbksb(matrix& u, vector<double>& w, matrix& v, vector<double>& b, vector<double>& x);
void svdcmp(matrix& a, vector<double>& w, matrix& v);

//-----------------------------------------------------------------------------
vector<double> operator / (vector<double>& b, matrix& m)
{
	int n = (int)b.size();

	vector<double> x(b);
	vector<int> indx(n);
	matrix a(m);

	ludcmp(a, n, &indx[0]);
	lubksb(a, n, &indx[0], &x[0]);

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

//-----------------------------------------------------------------------------
void matrix::alloc(int nr, int nc)
{
	m_nr = nr;
	m_nc = nc;
	m_nsize = nr*nc;

	m_pd = new double [m_nsize];
	m_pr = new double*[nr];
	for (int i=0; i<nr; i++) m_pr[i] = m_pd + i*nc;
}

//-----------------------------------------------------------------------------
//! Constructor for matrix class. 
matrix::matrix(int nr, int nc)
{
	alloc(nr, nc);
}

//-----------------------------------------------------------------------------
//! matrix destructor
void matrix::clear()
{
	if (m_pr) delete [] m_pr;
	if (m_pd) delete [] m_pd;
	m_pd = 0;
	m_pr = 0;
	m_nr = m_nc = 0;
}

//-----------------------------------------------------------------------------
//! Copy constructor for matrix class. 
matrix::matrix(const matrix& m)
{
	alloc(m.m_nr, m.m_nc);
	for (int i=0; i<m_nsize; ++i) m_pd[i] = m.m_pd[i];
}

//-----------------------------------------------------------------------------
//! constructor
matrix::matrix(const mat3d& m)
{
	alloc(3, 3);
	m_pr[0][0] = m[0][0]; m_pr[0][1] = m[0][1]; m_pr[0][2] = m[0][2];
	m_pr[1][0] = m[1][0]; m_pr[1][1] = m[1][1]; m_pr[1][2] = m[1][2];
	m_pr[2][0] = m[2][0]; m_pr[2][1] = m[2][1]; m_pr[2][2] = m[2][2];
}

//-----------------------------------------------------------------------------
//! assignment operator
matrix& matrix::operator = (const mat3d& m)
{
	resize(3, 3);
	m_pr[0][0] = m[0][0]; m_pr[0][1] = m[0][1]; m_pr[0][2] = m[0][2];
	m_pr[1][0] = m[1][0]; m_pr[1][1] = m[1][1]; m_pr[1][2] = m[1][2];
	m_pr[2][0] = m[2][0]; m_pr[2][1] = m[2][1]; m_pr[2][2] = m[2][2];
	return *this;
}

//-----------------------------------------------------------------------------
matrix& matrix::operator = (const matrix& m)
{
	if ((m.m_nr != m_nr) || (m.m_nc != m_nc))
	{
		clear();
		alloc(m.m_nr, m.m_nc);
	}
	for (int i=0; i<m_nsize; ++i) m_pd[i] = m.m_pd[i];

	return (*this);
}

//-----------------------------------------------------------------------------
void matrix::resize(int nr, int nc)
{
	if ((nr != m_nr) || (nc != m_nc))
	{
		clear();
		alloc(nr, nc);
	}
}

//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
void matrix::adds(int i, int j, const matrix& m, double s)
{
	int mr = m.rows();
	int mc = m.columns();
	for (int r = 0; r < mr; ++r)
		for (int c = 0; c < mc; ++c) 
			m_pr[i + r][j + c] += m[r][c]*s;
}

//-----------------------------------------------------------------------------
// Calculate the LU decomposition of this matrix. Note that this will modify
// the matrix. This is used for repeated solves of a linear system. Use 
// lusolve for solving after lufactor.
void matrix::lufactor(vector<int>& indx)
{
	// make sure this is a square matrix
	assert(m_nr == m_nc);

	// do a LU decomposition
	int n = m_nr;
	indx.resize(n);
	ludcmp(*(this), n, &indx[0]);
}

//-----------------------------------------------------------------------------
// Solve the linear system Ax=b, where A has been factored using lufactor.
// The indx array is the same one that was returned from lufactor
void matrix::lusolve(vector<double>& b, vector<int>& indx)
{
	// make sure this is a square matrix
	assert(m_nr == m_nc);
	lubksb(*(this), m_nr, &indx[0], &b[0]);
}

//-----------------------------------------------------------------------------
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
	ludcmp(a, n, &indx[0]);

	// allocate the inverse matrix
	matrix ai(n, n);

	// do a backsubstituation on the columns of a
	vector<double> b; b.assign(n, 0);
	for (int j=0; j<n; ++j)
	{
		b[j] = 1;
		lubksb(a, n, &indx[0], &b[0]);

		for (int i=0; i<n; ++i)
		{
			ai[i][j] = b[i];
			b[i] = 0;
		}
	}

	return ai;
}

//-----------------------------------------------------------------------------
// Matrix using singular value decomposition
matrix matrix::svd_inverse()
{
	matrix U(*this);
	matrix V(m_nr, m_nc);
	vector<double> w(m_nc);

	// calculate the decomposition
	svdcmp(U, w, V);

	matrix Ai(m_nc, m_nr); // inverse
	for (int i=0; i<m_nc; ++i)
		for (int j=0; j<m_nr; ++j)
		{
			double s = 0.0;
			for (int k=0;k<m_nc; ++k)
			{
				if (w[k] > 0.0)
				{
					s += (V[i][k]*U[j][k])/w[k];
				}
			}
			Ai[i][j] = s;
		}
	return Ai;
}

//-----------------------------------------------------------------------------
matrix matrix::transpose()
{
	int i, j;
	matrix At(m_nc, m_nr);
	for (i=0; i<m_nr; ++i)
		for (j=0; j<m_nc; ++j) At[j][i] = m_pr[i][j];
	return At;
}

//-----------------------------------------------------------------------------
matrix& matrix::operator += (const matrix& m)
{
	assert((m_nr == m.m_nr ) && (m_nc == m.m_nc));
	for (int i=0; i<m_nsize; ++i) m_pd[i] += m.m_pd[i];
	return (*this);
}

//-----------------------------------------------------------------------------
matrix& matrix::operator -= (const matrix& m)
{
	assert((m_nr == m.m_nr ) && (m_nc == m.m_nc));
	for (int i=0; i<m_nsize; ++i) m_pd[i] -= m.m_pd[i];
	return (*this);
}

//-----------------------------------------------------------------------------
// calculate outer product of a vector to produce a matrix
matrix outer_product(vector<double>& a)
{
	int n = (int) a.size();
	matrix m(n,n);
	for (int i=0; i<n; ++i)
	{
		for (int j=0; j<n; ++j) m[i][j] = a[i]*a[j];
	}
	return m;
}

//-----------------------------------------------------------------------------
// Calculates the infinity norm. That is, the max of the absolute row sum.
double matrix::inf_norm()
{
	matrix& self = (*this);
	double m = 0;
	for (int j=0; j<m_nr; ++j)
	{
		double s = 0;
		for (int i=0; i<m_nc; ++i) s += fabs(self[j][i]);
		if (s > m) m = s;
	}
	return m;
}

//-----------------------------------------------------------------------------
void matrix::get(int i, int j, int rows, int cols, matrix& A) const
{
	// make sure we create a valid matrix
	if ((rows <= 0) || (cols <= 0)) return;

	// initialize the matrix
	A.resize(rows, cols);
	A.zero();

	// make sure the bounds are within this matrix
	if ((i >= m_nr) || (j >= m_nc)) return;
	if ((i + rows <= 0) || (j + cols <= 0)) return;

	// set range
	int r0 = 0;
	int r1 = rows - 1;
	int c0 = 0;
	int c1 = cols - 1;
	if (i < 0) r0 -= i;
	if (j < 0) c0 -= j;
	if (i + rows > m_nr) r1 -= rows + i - m_nr;
	if (j + cols > m_nc) c1 -= cols + j - m_nc;

	for (int r=r0; r<=r1; ++r)
		for (int c=c0; c<=c1; ++c)
				A[r][c] = (*this)(i+r, j+c);
}

//-----------------------------------------------------------------------------
// fill a matrix
void matrix::fill(int i, int j, int rows, int cols, double val)
{
	if ((i >= m_nr) || (j >= m_nc)) return;
	if ((i + rows <= 0) || (j + cols <= 0)) return;

	int r0 = 0;
	int r1 = rows - 1;
	int c0 = 0;
	int c1 = cols - 1;
	if (i < 0) r0 -= i;
	if (j < 0) c0 -= j;
	if (i + rows > m_nr) r1 -= rows + i - m_nr;
	if (j + cols > m_nc) c1 -= cols + j - m_nc;

	for (int r = r0; r<=r1; ++r)
		for (int c = c0; c<=c1; ++c)
			(*this)(i + r, j + c) = val;
}

//-----------------------------------------------------------------------------
// solve the linear system Ax=b
void matrix::solve(const vector<double>& b, vector<double>& x)
{
	matrix A(*this);
	vector<int> index;
	A.lufactor(index);
	x = b;
	A.lusolve(x, index);
}

//-----------------------------------------------------------------------------
matrix matrix::operator + (const matrix& m)
{
	matrix s(*this);
	for (int i=0; i<m_nr; ++i)
		for (int j=0; j<m_nc; ++j)
			s(i,j) += m(i,j);

	return s;
}

//-----------------------------------------------------------------------------
matrix matrix::operator - (const matrix& m)
{
	matrix s(*this);
	for (int i = 0; i<m_nr; ++i)
		for (int j = 0; j<m_nc; ++j)
			s(i, j) -= m(i, j);

	return s;
}
