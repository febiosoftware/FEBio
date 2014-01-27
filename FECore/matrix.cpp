// matrix.cpp: implementation of the matrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "matrix.h"
#include <assert.h>

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
	int n = b.size();

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
//! Constructor for matrix class. 
matrix::matrix(int nr, int nc)
{
	m_nr = nr;
	m_nc = nc;
	m_nsize = nr*nc;

	m_pd = new double [m_nsize];
	m_pr = new double*[nr];

	for (int i=0; i<nr; i++) m_pr[i] = m_pd + i*nc;
}

//-----------------------------------------------------------------------------
//! matrix destructor
matrix::~matrix()
{
	delete [] m_pd;
	delete [] m_pr;
}

//-----------------------------------------------------------------------------
//! Copy constructor for matrix class. 
matrix::matrix(const matrix& m)
{
	m_nr = m.m_nr;
	m_nc = m.m_nc;
	m_nsize = m_nr*m_nc;

	m_pd = new double[m_nsize];
	m_pr = new double*[m_nr];

	int i;
	for (i=0; i<m_nr; ++i) m_pr[i] = m_pd + i*m_nc;
	for (i=0; i<m_nsize; ++i) m_pd[i] = m.m_pd[i];
}

//-----------------------------------------------------------------------------
matrix& matrix::operator = (const matrix& m)
{
	int i;
	if ((m.m_nr != m_nr) || (m.m_nc != m_nc))
	{
		delete [] m_pd;
		delete [] m_pr;

		m_nr = m.m_nr;
		m_nc = m.m_nc;
		m_nsize = m_nr*m_nc;

		m_pd = new double[m_nsize];
		m_pr = new double*[m_nr];

		for (i=0; i<m_nr; ++i) m_pr[i] = m_pd + i*m_nc;
	}

	for (i=0; i<m_nsize; ++i) m_pd[i] = m.m_pd[i];

	return (*this);
}

//-----------------------------------------------------------------------------
void matrix::resize(int nr, int nc)
{
	if ((nr != m_nr) || (nc != m_nc))
	{
		m_nr = nr;
		m_nc = nc;
		m_nsize = nr*nc;

		if (m_pd) delete [] m_pd;
		if (m_pr) delete [] m_pr;

		m_pd = new double [m_nsize];
		m_pr = new double*[nr];

		for (int i=0; i<nr; i++) m_pr[i] = m_pd + i*nc;
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
