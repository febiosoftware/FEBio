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
#include "mat3d.h"
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
	~matrix() { clear(); }

	//! access operator
	double * operator [] (int l) { return m_pr[l]; }
	double& operator () (int i, int j) { return m_pr[i][j]; }
	double operator () (int i, int j) const { return m_pr[i][j]; }
	operator double** () { return m_pr; }

	int rows   () const { return m_nr; }
	int columns() const { return m_nc; }

	void zero() { memset(m_pd, 0, sizeof(double)*m_nsize); }

	//! matrix transpose
	matrix transpose();

	//! matrix inversion
	matrix inverse();

	//! matrix inverse using SVD
	matrix svd_inverse();

	//! matrix operators
	matrix operator * (const matrix& m);

	matrix& operator += (const matrix& m);

	matrix& operator -= (const matrix& m);

	matrix& operator *=(double g)
	{
		for (int i=0; i<m_nsize; ++i) m_pd[i] *= g;
		return *this;
	}

	// calculate the LU decomposition
	// note that this modifies the matrix
	void lufactor(vector<int>& indx);

	// solve using the lu factor calculated with lufactor
	void lusolve(vector<double>& b, vector<int>& indx);

	// infinity-norm
	double inf_norm();

public:
	void set(int i, int j, const mat3d& a)
	{
		m_pr[i][j] = a(0,0); m_pr[i][j+1] = a(0,1); m_pr[i][j+2] = a(0,2); i++;
		m_pr[i][j] = a(1,0); m_pr[i][j+1] = a(1,1); m_pr[i][j+2] = a(1,2); i++;
		m_pr[i][j] = a(2,0); m_pr[i][j+1] = a(2,1); m_pr[i][j+2] = a(2,2);
	}

	void add(int i, int j, const mat3d& a)
	{
		m_pr[i][j] += a(0,0); m_pr[i][j+1] += a(0,1); m_pr[i][j+2] += a(0,2); i++;
		m_pr[i][j] += a(1,0); m_pr[i][j+1] += a(1,1); m_pr[i][j+2] += a(1,2); i++;
		m_pr[i][j] += a(2,0); m_pr[i][j+1] += a(2,1); m_pr[i][j+2] += a(2,2);
	}

	void sub(int i, int j, const mat3d& a)
	{
		m_pr[i][j] -= a(0,0); m_pr[i][j+1] -= a(0,1); m_pr[i][j+2] -= a(0,2); i++;
		m_pr[i][j] -= a(1,0); m_pr[i][j+1] -= a(1,1); m_pr[i][j+2] -= a(1,2); i++;
		m_pr[i][j] -= a(2,0); m_pr[i][j+1] -= a(2,1); m_pr[i][j+2] -= a(2,2);
	}

	// copy-lower-triangular
	// make the matrix symmetric by copying the lower triangular part
	void copy_lt()
	{
		assert(m_nr==m_nc);
		if (m_nr != m_nc) return;
		for (int i=0; i<m_nr; ++i)
		{
			for (int j=i+1; j<m_nr; ++j)
			{
				m_pr[i][j] = m_pr[j][i];
			}
		}
	}

	// copy-upper-triangular
	// make the matrix symmetric by copying the upper triangular part
	void copy_ut()
	{
		assert(m_nr==m_nc);
		if (m_nr != m_nc) return;
		for (int i=0; i<m_nr; ++i)
		{
			for (int j=i+1; j<m_nr; ++j)
			{
				m_pr[j][i] = m_pr[i][j];
			}
		}
	}


private:
	void alloc(int nr, int nc);
	void clear();

protected:
	double**	m_pr;	// pointer to rows
	double*		m_pd;	// matrix elements

	int	m_nr;		// nr of rows
	int	m_nc;		// nr of columns
	int	m_nsize;	// size of matrix (ie. total nr of elements = nr*nc)
};

vector<double> operator / (vector<double>& b, matrix& m);
vector<double> operator * (matrix& m, vector<double>& b);
matrix outer_product(vector<double>& a);

#endif // !defined(AFX_MATRIX_H__C0F2C6F6_AE26_4C7F_8C70_5A7BF5DD421E__INCLUDED_)
