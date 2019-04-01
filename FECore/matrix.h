/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

#pragma once
#include <memory.h>
#include <vector>
#include "fecore_api.h"
#include "mat3d.h"
using namespace std;

//-----------------------------------------------------------------------------
//! General purpose matrix class.
class FECORE_API matrix
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
	const double* operator [] (int l) const { return m_pr[l]; }
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

	matrix operator + (const matrix& m);

	matrix operator - (const matrix& m);

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

	// solve the linear system Ax=b
	void solve(const vector<double>& b, vector<double>& x);

	// infinity-norm
	double inf_norm();

public:
	void set(int i, int j, const mat3d& a);
	
	void add(int i, int j, const mat3ds& a);
	void add(int i, int j, const mat3da& a);
	void add(int i, int j, const mat3dd& a);
	void add(int i, int j, const mat3d&  a);

	void sub(int i, int j, const mat3ds& a);
	void sub(int i, int j, const mat3da& a);
	void sub(int i, int j, const mat3dd& a);
	void sub(int i, int j, const mat3d&  a);

	void get(int i, int j, mat3d& a);

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

	// extract a matrix block
	// the returned matrix will have the dimensions rows x cols
	// if the matrix doesn't fit in this matrix, the missing entries will be set to zero
	void get(int i, int j, int rows, int cols, matrix& A) const;

	// fill a matrix
	void fill(int i, int j, int rows, int cols, double val);

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

vector<double> FECORE_API operator / (vector<double>& b, matrix& m);
vector<double> FECORE_API operator * (matrix& m, vector<double>& b);
matrix FECORE_API outer_product(vector<double>& a);

inline void matrix::set(int i, int j, const mat3d& a)
{
	m_pr[i][j] = a(0,0); m_pr[i][j+1] = a(0,1); m_pr[i][j+2] = a(0,2); i++;
	m_pr[i][j] = a(1,0); m_pr[i][j+1] = a(1,1); m_pr[i][j+2] = a(1,2); i++;
	m_pr[i][j] = a(2,0); m_pr[i][j+1] = a(2,1); m_pr[i][j+2] = a(2,2);
}

inline void matrix::add(int i, int j, const mat3ds& a)
{
	m_pr[i][j] += a.xx(); m_pr[i][j+1] += a.xy(); m_pr[i][j+2] += a.xz(); i++;
	m_pr[i][j] += a.xy(); m_pr[i][j+1] += a.yy(); m_pr[i][j+2] += a.yz(); i++;
	m_pr[i][j] += a.xz(); m_pr[i][j+1] += a.yz(); m_pr[i][j+2] += a.zz();
}

inline void matrix::add(int i, int j, const mat3da& a)
{
	m_pr[i][j+1] += a.xy(); m_pr[i][j+2] += a.xz(); i++;
	m_pr[i][j  ] -= a.xy(); m_pr[i][j+2] += a.yz(); i++;
	m_pr[i][j  ] -= a.xz(); m_pr[i][j+1] -= a.yz();
}

inline void matrix::add(int i, int j, const mat3dd& a)
{
	m_pr[i][j  ] += a.diag(0); i++;
	m_pr[i][j+1] += a.diag(1); i++;
	m_pr[i][j+2] += a.diag(2);
}

inline void matrix::add(int i, int j, const mat3d& a)
{
	m_pr[i][j] += a(0,0); m_pr[i][j+1] += a(0,1); m_pr[i][j+2] += a(0,2); i++;
	m_pr[i][j] += a(1,0); m_pr[i][j+1] += a(1,1); m_pr[i][j+2] += a(1,2); i++;
	m_pr[i][j] += a(2,0); m_pr[i][j+1] += a(2,1); m_pr[i][j+2] += a(2,2);
}

inline void matrix::sub(int i, int j, const mat3ds& a)
{
	m_pr[i][j] -= a.xx(); m_pr[i][j+1] -= a.xy(); m_pr[i][j+2] -= a.xz(); i++;
	m_pr[i][j] -= a.xy(); m_pr[i][j+1] -= a.yy(); m_pr[i][j+2] -= a.yz(); i++;
	m_pr[i][j] -= a.xz(); m_pr[i][j+1] -= a.yz(); m_pr[i][j+2] -= a.zz();
}

inline void matrix::sub(int i, int j, const mat3da& a)
{
	m_pr[i][j+1] -= a.xy(); m_pr[i][j+2] -= a.xz(); i++;
	m_pr[i][j  ] += a.xy(); m_pr[i][j+2] -= a.yz(); i++;
	m_pr[i][j  ] += a.xz(); m_pr[i][j+1] += a.yz();
}

inline void matrix::sub(int i, int j, const mat3dd& a)
{
	m_pr[i][j  ] -= a.diag(0); i++;
	m_pr[i][j+1] -= a.diag(1); i++;
	m_pr[i][j+2] -= a.diag(2);
}

inline void matrix::sub(int i, int j, const mat3d& a)
{
	m_pr[i][j] -= a(0,0); m_pr[i][j+1] -= a(0,1); m_pr[i][j+2] -= a(0,2); i++;
	m_pr[i][j] -= a(1,0); m_pr[i][j+1] -= a(1,1); m_pr[i][j+2] -= a(1,2); i++;
	m_pr[i][j] -= a(2,0); m_pr[i][j+1] -= a(2,1); m_pr[i][j+2] -= a(2,2);
}

inline void matrix::get(int i, int j, mat3d& a)
{
	a[0][0] = m_pr[i  ][j]; a[0][1] = m_pr[i  ][j+1]; a[0][2] = m_pr[i  ][j+2];
	a[1][0] = m_pr[i+1][j]; a[1][1] = m_pr[i+1][j+1]; a[1][2] = m_pr[i+1][j+2];
	a[2][0] = m_pr[i+2][j]; a[2][1] = m_pr[i+2][j+1]; a[2][2] = m_pr[i+2][j+2];
}
