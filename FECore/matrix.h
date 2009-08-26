// matrix.h: interface for the matrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIX_H__C0F2C6F6_AE26_4C7F_8C70_5A7BF5DD421E__INCLUDED_)
#define AFX_MATRIX_H__C0F2C6F6_AE26_4C7F_8C70_5A7BF5DD421E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <memory.h>
#include "vector.h"

class matrix  
{
public:
	matrix() : m_nr(0), m_nc(0), m_nsize(0), m_pd(0), m_pr(0) {}

	matrix(int nr, int nc)
	{
		m_nr = nr;
		m_nc = nc;
		m_nsize = nr*nc;

		m_pd = new double [m_nsize];
		m_pr = new double*[nr];

		for (int i=0; i<nr; i++) m_pr[i] = m_pd + i*nc;
	}
	matrix(const matrix& m)
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

	matrix& operator = (const matrix& m)
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

	void Create(int nr, int nc)
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

	~matrix()
	{
		delete [] m_pd;
		delete [] m_pr;
	}

	double * operator [] (int l) { return m_pr[l]; }

	double& operator () (int i, int j) { return m_pr[i][j]; }

	double operator () (int i, int j) const { return m_pr[i][j]; }

	int rows   () { return m_nr; }
	int columns() { return m_nc; }

	void zero() { memset(m_pd, 0, sizeof(double)*m_nsize); }

	operator double** () { return m_pr; }

	matrix transpose()
	{
		int i, j;
		matrix At(m_nc, m_nr);
		for (i=0; i<m_nr; ++i)
			for (j=0; j<m_nc; ++j) At[j][i] = m_pr[i][j];
		return At;
	}

	matrix inverse();

	matrix operator * (const matrix& m);

protected:
	double**	m_pr;	// pointer to rows
	double*		m_pd;	// matrix elements

	int	m_nr;		// nr of rows
	int	m_nc;		// nr of columns
	int	m_nsize;	// size of matrix (ie. total nr of elements = nr*nc)
};

vector<double> operator / (vector<double>& b, matrix& m);

#endif // !defined(AFX_MATRIX_H__C0F2C6F6_AE26_4C7F_8C70_5A7BF5DD421E__INCLUDED_)
