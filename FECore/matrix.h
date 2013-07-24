// matrix.h: interface for the matrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIX_H__C0F2C6F6_AE26_4C7F_8C70_5A7BF5DD421E__INCLUDED_)
#define AFX_MATRIX_H__C0F2C6F6_AE26_4C7F_8C70_5A7BF5DD421E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <memory.h>
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! General purpose matrix class.
class matrix  
{
public:
	//! constructor
	matrix() : m_nr(0), m_nc(0), m_nsize(0), m_pd(0), m_pr(0) {}

	//! constructor
	matrix(int nr, int nc);

	//! copy constructor
	matrix(const matrix& m);

	//! assignment operator
	matrix& operator = (const matrix& m);

	//! Matrix reallocation
	void resize(int nr, int nc);

	//! destructor
	~matrix();

	//! access operator
	double * operator [] (int l) { return m_pr[l]; }
	double& operator () (int i, int j) { return m_pr[i][j]; }
	double operator () (int i, int j) const { return m_pr[i][j]; }
	operator double** () { return m_pr; }

	int rows   () { return m_nr; }
	int columns() { return m_nc; }

	void zero() { memset(m_pd, 0, sizeof(double)*m_nsize); }

	//! matrix transpose
	matrix transpose();

	//! matrix inversion
	matrix inverse();

	//! matrix operators
	matrix operator * (const matrix& m);

	matrix& operator += (const matrix& m);

	matrix& operator -= (const matrix& m);

protected:
	double**	m_pr;	// pointer to rows
	double*		m_pd;	// matrix elements

	int	m_nr;		// nr of rows
	int	m_nc;		// nr of columns
	int	m_nsize;	// size of matrix (ie. total nr of elements = nr*nc)
};

vector<double> operator / (vector<double>& b, matrix& m);
vector<double> operator * (matrix& m, vector<double>& b);

#endif // !defined(AFX_MATRIX_H__C0F2C6F6_AE26_4C7F_8C70_5A7BF5DD421E__INCLUDED_)
