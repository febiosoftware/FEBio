#include "stdafx.h"
#include "CSRMatrix.h"
#include <assert.h>
#include "mkl_spblas.h"

CSRMatrix::CSRMatrix() : m_nr(0), m_nc(0)
{
}

// create a matrix of given size
CSRMatrix::CSRMatrix(int rows, int cols) : m_nr(rows), m_nc(cols)
{
	m_rowIndex.resize(rows+1, 0);
}

// Create matrix
void CSRMatrix::create(int nr, int nc)
{
	m_nr = nr;
	m_nc = nc;
	m_rowIndex.resize(nr+1, 0);
	m_columns.clear();
	m_values.clear();
}

// copy constructor
CSRMatrix::CSRMatrix(const CSRMatrix& A)
{
	m_nr = A.m_nr;
	m_nc = A.m_nc;
	m_rowIndex = A.m_rowIndex;
	m_columns = A.m_columns;
	m_values = A.m_values;
}

// assignment operator
void CSRMatrix::operator = (const CSRMatrix& A)
{
	m_nr = A.m_nr;
	m_nc = A.m_nc;
	m_rowIndex = A.m_rowIndex;
	m_columns = A.m_columns;
	m_values = A.m_values;
}

// sets a value, inserting it if necessary
void CSRMatrix::set(int i, int j, double val)
{
	assert((i >= 0) && (i < m_nr));
	assert((j >= 0) && (j < m_nc));

	// get the start column index for the given row and the non-zero count for that row
	int col = m_rowIndex[i];
	int count = m_rowIndex[i+1] - m_rowIndex[i];

	// see if the row is empty
	if (count == 0)
	{
		m_columns.insert(m_columns.begin() + col, j);
		m_values.insert(m_values.begin() + col, val);
	}
	else
	{
		// see if this column would be before the first entry
		if (j < m_columns[col])
		{
			m_columns.insert(m_columns.begin() + col, j);
			m_values.insert(m_values.begin() + col, val);
		}
		// see if this column would be past the last entry
		else if (j > m_columns[col + count - 1])
		{
			m_columns.insert(m_columns.begin() + col + count, j);
			m_values.insert(m_values.begin() + col + count, val);
		}
		else {
			// find the column index
			for (int n=0; n<count; ++n)
			{
				// see if it alreay exists
				if (m_columns[col+n] == j)
				{
					// if so, replace the value and return
					m_values[col+n] = val;
					return;
				}
				else if (m_columns[col + n] > j)
				{
					// we found an index that is bigger so insert this value before this one
					m_columns.insert(m_columns.begin() + col + n, j);
					m_values.insert(m_values.begin() + col + n, val);
					break;
				}
			}
		}
	}

	// increase row counts
	for (int n = i + 1; n <= m_nr; ++n) m_rowIndex[n]++;
	assert(m_rowIndex[m_nr] == m_values.size());
	assert(m_columns.size() == m_values.size());
}

// get a value
double CSRMatrix::operator () (int i, int j) const
{
	assert((i >= 0) && (i < m_nr));
	assert((j >= 0) && (j < m_nc));

	int col = m_rowIndex[i];
	int count = m_rowIndex[i + 1] - m_rowIndex[i];

	if (count == 0) return 0.0;
	if (j < m_columns[col]) return 0.0;
	if (j > m_columns[col + count - 1]) return 0.0;
	for (int n=0; n<count; ++n)
	{
		if (m_columns[col+n] == j) return m_values[col+n];
	}
	return 0.0;
}

// see if a matrix entry was allocated
bool CSRMatrix::isAlloc(int i, int j) const
{
	assert((i >= 0) && (i < m_nr));
	assert((j >= 0) && (j < m_nc));

	int col = m_rowIndex[i];
	int count = m_rowIndex[i + 1] - m_rowIndex[i];

	if (count == 0) return false;
	if (j < m_columns[col]) return false;
	if (j > m_columns[col + count - 1]) return false;
	for (int n = 0; n<count; ++n)
	{
		if (m_columns[col + n] == j) return true;
	}
	return false;
}

CSRMatrix CSRMatrix::operator + (const CSRMatrix& A)
{
	assert((A.m_nr == m_nr) && (A.m_nc == m_nc));

	// create empty matrix
	CSRMatrix S(m_nr, m_nc);

	// count the number of non-zero values we'll have
	int nnz = 0;
	for (int i=0; i<m_nr; ++i)
	{
		int c0 = m_rowIndex[i];
		int c1 = A.m_rowIndex[i];

		int cn0 = m_rowIndex[i+1] - c0;
		int cn1 = A.m_rowIndex[i+1] - c1;

		int n0 = 0;
		int n1 = 0;
		while ((n0 < cn0) && (n1 < cn1))
		{
			int m0 = m_columns[c0 + n0];
			int m1 = A.m_columns[c1 + n1];
			if (m0 < m1) n0++;
			else if (m1 < m0) n1++;
			else if (m0 == m1)
			{
				n0++;
				n1++;
			}
			nnz++;
		}
		nnz += cn0 - n0;
		nnz += cn1 - n1;

		S.m_rowIndex[i+1] = nnz;
	}

	// allocate data
	S.m_values.resize(nnz);
	S.m_columns.resize(nnz);

	// fill in the values
	nnz = 0;
	for (int i = 0; i<m_nr; ++i)
	{
		int c0 = m_rowIndex[i];
		int c1 = A.m_rowIndex[i];

		int cn0 = m_rowIndex[i + 1] - c0;
		int cn1 = A.m_rowIndex[i + 1] - c1;

		int n0 = 0;
		int n1 = 0;
		while ((n0 < cn0) && (n1 < cn1))
		{
			int m0 = m_columns[c0 + n0];
			int m1 = A.m_columns[c1 + n1];

			double v0 = m_values[c0 + n0];
			double v1 = A.m_values[c1 + n1];

			if (m0 < m1) 
			{
				S.m_columns[nnz] = m0;
				S.m_values[nnz] = v0;
				n0++;
			}
			else if (m1 < m0) 
			{
				S.m_columns[nnz] = m1;
				S.m_values[nnz] = v1;
				n1++;
			}
			else if (m0 == m1)
			{
				S.m_columns[nnz] = m0;
				S.m_values[nnz] = v0 + v1;
				n0++;
				n1++;
			}
			nnz++;
		}

		if (n0 < cn0)
		{
			assert(n1 == cn1);
			for (int n=n0; n<cn0; ++n)
			{
				S.m_columns[nnz] = m_columns[c0 + n];
				S.m_values[nnz] = m_values[c0 + n];
				nnz++;
			}
		}
		if (n1 < cn1)
		{
			assert(n0 == cn0);
			for (int n = n1; n<cn1; ++n)
			{
				S.m_columns[nnz] = A.m_columns[c1 + n];
				S.m_values[nnz] = A.m_values[c1 + n];
				nnz++;
			}
		}
	}

	return S;
}

CSRMatrix CSRMatrix::operator - (const CSRMatrix& A)
{
	assert((A.m_nr == m_nr) && (A.m_nc == m_nc));

	// create empty matrix
	CSRMatrix S(m_nr, m_nc);

	// count the number of non-zero values we'll have
	int nnz = 0;
	for (int i = 0; i<m_nr; ++i)
	{
		int c0 = m_rowIndex[i];
		int c1 = A.m_rowIndex[i];

		int cn0 = m_rowIndex[i + 1] - c0;
		int cn1 = A.m_rowIndex[i + 1] - c1;

		int n0 = 0;
		int n1 = 0;
		while ((n0 < cn0) && (n1 < cn1))
		{
			int m0 = m_columns[c0 + n0];
			int m1 = A.m_columns[c1 + n1];
			if (m0 < m1) n0++;
			else if (m1 < m0) n1++;
			else if (m0 == m1)
			{
				n0++;
				n1++;
			}
			nnz++;
		}
		nnz += cn0 - n0;
		nnz += cn1 - n1;

		S.m_rowIndex[i + 1] = nnz;
	}

	// allocate data
	S.m_values.resize(nnz);
	S.m_columns.resize(nnz);

	// fill in the values
	nnz = 0;
	for (int i = 0; i<m_nr; ++i)
	{
		int c0 = m_rowIndex[i];
		int c1 = A.m_rowIndex[i];

		int cn0 = m_rowIndex[i + 1] - c0;
		int cn1 = A.m_rowIndex[i + 1] - c1;

		int n0 = 0;
		int n1 = 0;
		while ((n0 < cn0) && (n1 < cn1))
		{
			int m0 = m_columns[c0 + n0];
			int m1 = A.m_columns[c1 + n1];

			double v0 = m_values[c0 + n0];
			double v1 = A.m_values[c1 + n1];

			if (m0 < m1)
			{
				S.m_columns[nnz] = m0;
				S.m_values[nnz] = v0;
				n0++;
			}
			else if (m1 < m0)
			{
				S.m_columns[nnz] = m1;
				S.m_values[nnz] = -v1;
				n1++;
			}
			else if (m0 == m1)
			{
				S.m_columns[nnz] = m0;
				S.m_values[nnz] = v0 - v1;
				n0++;
				n1++;
			}
			nnz++;
		}

		if (n0 < cn0)
		{
			assert(n1 == cn1);
			for (int n = n0; n<cn0; ++n)
			{
				S.m_columns[nnz] = m_columns[c0 + n];
				S.m_values[nnz] = m_values[c0 + n];
				nnz++;
			}
		}
		if (n1 < cn1)
		{
			assert(n0 == cn0);
			for (int n = n1; n<cn1; ++n)
			{
				S.m_columns[nnz] = A.m_columns[c1 + n];
				S.m_values[nnz] = -A.m_values[c1 + n];
				nnz++;
			}
		}
	}

	return S;
}

std::vector<double> CSRMatrix::operator * (const std::vector<double>& a)
{
	assert(m_nr = (int)a.size());
	std::vector<double> b(m_nr);
	for (int i=0; i<m_nr; ++i)
	{
		int col = m_rowIndex[i];
		int cn = m_rowIndex[i+1] - col;

		double s = 0.0;
		for (int j=0; j<cn; ++j)
		{
			s += m_values[col + j]*a[ m_columns[col + j] ];
		}
		b[i] = s;
	}
	return b;
}

void CSRMatrix::multv(const std::vector<double>& x, std::vector<double>& r)
{
/*	for (int i = 0; i<m_nr; ++i)
	{
		int col = m_rowIndex[i];
		int cn = m_rowIndex[i + 1] - col;

		double s = 0.0;
		for (int j = 0; j<cn; ++j)
		{
			s += m_values[col + j] * x[m_columns[col + j]];
		}
		r[i] = s;
	}
*/
	mkl_cspblas_dcsrgemv("N", &m_nr, &m_values[0], &m_rowIndex[0], &m_columns[0], (double*) &x[0], &r[0]);

}

void CSRMatrix::multv(const double* x, double* r)
{
/*	for (int i = 0; i<m_nr; ++i)
	{
		int col = m_rowIndex[i];
		int cn = m_rowIndex[i + 1] - col;

		double s = 0.0;
		for (int j = 0; j<cn; ++j)
		{
			s += m_values[col + j] * x[m_columns[col + j]];
		}
		r[i] = s;
	}
*/
	mkl_cspblas_dcsrgemv("N", &m_nr, &m_values[0], &m_rowIndex[0], &m_columns[0], (double*)x, r);
}

// normalize the matrix 
void CSRMatrix::normalize(const std::vector<double>& l, const std::vector<double>& r)
{
	for (int i=0; i<m_nr; ++i)
	{
		int col = m_rowIndex[i];
		int cn = m_rowIndex[i + 1] - col;

		for (int j=0; j<cn; ++j)
		{
			m_values[col + j] *= l[i]*r[ m_columns[col + j ]];
		}
	}
}
