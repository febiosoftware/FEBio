/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "MMatrix.h"
using namespace std;

//-----------------------------------------------------------------------------
MMatrix::MMatrix() : MItem(MMATRIX)
{
	m_nrow = 0; m_ncol = 0;
}

//-----------------------------------------------------------------------------
MMatrix::~MMatrix()
{
	for (int i=0;i<m_nrow; ++i)
		for (int j=0; j<m_ncol; ++j) delete m_d[i][j];
	m_d.clear();
	m_nrow = 0;
	m_ncol = 0;
}

//-----------------------------------------------------------------------------
void MMatrix::Create(int nrow, int ncol)
{
	m_d.clear();
	m_d.assign(nrow, vector<MItem*>(ncol));
	m_nrow = nrow;
	m_ncol = ncol;
}

//-----------------------------------------------------------------------------
MItem* MMatrix::copy() const
{
	MMatrix* pm = new MMatrix();
	MMatrix& m = *pm;
	const MMatrix& a = (*this);
	m.Create(m_nrow, m_ncol);
	for (int i=0; i<m_nrow; ++i)
		for (int j=0; j<m_ncol; ++j) m[i][j] = a[i][j]->copy();
	return pm;
}

//-----------------------------------------------------------------------------
MMatrix* operator + (const MMatrix& A, const MMatrix& B)
{
	if ((A.rows() != B.rows()) || (A.columns() != B.columns())) throw InvalidOperation();
	MMatrix& C = *(new MMatrix());
	int N = A.rows();
	int M = A.columns();
	C.Create(A.rows(), A.columns());
	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<M; ++j)
		{
			MITEM aij = MITEM(A[i][j]->copy());
			MITEM bij = MITEM(B[i][j]->copy());
			C[i][j] = (aij + bij).copy();
		}
	}
	return &C;
}

//-----------------------------------------------------------------------------
MMatrix* operator - (const MMatrix& A, const MMatrix& B)
{
	if ((A.rows() != B.rows()) || (A.columns() != B.columns())) throw InvalidOperation();
	MMatrix& C = *(new MMatrix());
	int N = A.rows();
	int M = A.columns();
	C.Create(A.rows(), A.columns());
	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<M; ++j)
		{
			MITEM aij = MITEM(A[i][j]->copy());
			MITEM bij = MITEM(B[i][j]->copy());
			C[i][j] = (aij - bij).copy();
		}
	}
	return &C;
}

//-----------------------------------------------------------------------------
MMatrix* operator * (const MMatrix& A, const MMatrix& B)
{
	int RA = A.rows();
	int CA = A.columns();
	if (B.rows() != CA) throw InvalidOperation();
	int CB = B.columns();
	MMatrix& C = *(new MMatrix());
	C.Create(RA, CB);
	for (int i=0; i<RA; ++i)
	{
		for (int j=0; j<CB; ++j)
		{
			MITEM ai0 = A[i][0]->copy();
			MITEM b0j = B[0][j]->copy();
			MITEM cij = ai0*b0j;
			for (int k=1; k<CA; ++k)
			{
				MITEM aik = A[i][k]->copy();
				MITEM bkj = B[k][j]->copy();
				cij = cij + aik*bkj;
			}
			C[i][j] = cij.copy();
		}
	}
	return (&C);
}

//-----------------------------------------------------------------------------
MMatrix* operator * (const MMatrix& A, const MITEM& n)
{
	if (n == 1.0) return static_cast<MMatrix*>(A.copy());
	MMatrix& C = *(new MMatrix());
	int N = A.rows();
	int M = A.columns();
	C.Create(A.rows(), A.columns());
	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<M; ++j)
		{
			MITEM aij(A[i][j]->copy());
			MITEM f(n.copy());
			C[i][j] = (f*aij).copy();
		}
	}
	return &C;
}

//-----------------------------------------------------------------------------
MMatrix* operator / (const MMatrix& A, const MITEM& n)
{
	if (n == 1.0) return static_cast<MMatrix*>(A.copy());
	MMatrix& C = *(new MMatrix());
	int N = A.rows();
	int M = A.columns();
	C.Create(A.rows(), A.columns());
	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<M; ++j)
		{
			MITEM aij(A[i][j]->copy());
			MITEM f(n.copy());
			C[i][j] = (aij/f).copy();
		}
	}
	return &C;
}


//-----------------------------------------------------------------------------
MItem* matrix_transpose(const MMatrix& A)
{
	MMatrix& B = *(new MMatrix());
	B.Create(A.columns(), A.rows());
	for (int i=0; i<B.rows(); ++i)
		for (int j=0; j<B.columns(); ++j)
		{
			B[i][j] = A[j][i]->copy();
		}
	return &B;
}

//-----------------------------------------------------------------------------
MItem* matrix_trace(const MMatrix& A)
{
	if (A.rows() != A.columns()) throw InvalidOperation();
	MITEM t = A[0][0]->copy();
	for (int i=1; i<A.rows(); ++i)
	{
		MITEM aii = A[i][i]->copy();
		t = t + aii;
	}
	return t.copy();
}

//-----------------------------------------------------------------------------
MItem* matrix_determinant(const MMatrix& A)
{
	// make sure the matrix is square
	if (A.rows() != A.columns()) throw InvalidOperation();

	if (A.rows() == 1) return A(0,0)->copy();
	if (A.rows() == 2)
	{
		MITEM D = MITEM(A(0,0))*MITEM(A(1,1)) - MITEM(A(0,1))*MITEM(A(1,0));
		return D.copy();
	}
	if (A.rows() == 3)
	{
		MITEM D(0.0);
		D = D + MITEM(A(0,0))*MITEM(A(1,1))*MITEM(A(2,2));
		D = D + MITEM(A(0,1))*MITEM(A(1,2))*MITEM(A(2,0));
		D = D + MITEM(A(0,2))*MITEM(A(1,0))*MITEM(A(2,1));
		D = D - MITEM(A(0,2))*MITEM(A(1,1))*MITEM(A(2,0));
		D = D - MITEM(A(0,1))*MITEM(A(1,0))*MITEM(A(2,2));
		D = D - MITEM(A(0,0))*MITEM(A(1,2))*MITEM(A(2,1));
		return D.copy();
	}
	throw InvalidOperation();
	return 0;
}

//-----------------------------------------------------------------------------
MItem* matrix_contract(const MMatrix& A, const MMatrix& B)
{
	if ((A.rows() != B.rows())||(A.columns() != B.columns())) throw InvalidOperation();
	MITEM s(0.0);
	for (int i=0; i<A.rows(); ++i)
		for (int j=0; j<A.columns(); ++j)
		{
			MITEM aij = A[i][j]->copy();
			MITEM bij = B[i][j]->copy();
			s = s + aij*bij;
		}
	return s.copy();
}

//-----------------------------------------------------------------------------
MMatrix* matrix_identity(int n)
{
	MMatrix& m = *(new MMatrix());
	m.Create(n, n);
	for (int i=0; i<n; ++i)
		for (int j=0; j<n; ++j)
		{
			m[i][j] = new MConstant((i==j?1.0:0.0));
		}
	return &m;
}

//-----------------------------------------------------------------------------
//! Calculate the inverse of a matrix
MItem* matrix_inverse(const MMatrix& m)
{
	// make sure it is a square matrix
	int nrows = m.rows();
	int ncols = m.columns();
	if (nrows != ncols) throw InvalidOperation();

	// calculate the inverse, but we can only handle matrix up to order 3 for now
	if (nrows == 1)
	{
		MMatrix* pmi = new MMatrix;
		MMatrix& mi = *pmi;
		mi.Create(1,1);
		MITEM mij(m[0][0]->copy());
		mi[0][0] = (MITEM(1.0) / mij).copy();
		return pmi;
	}
	else if (nrows == 2)
	{
		MITEM m00(m[0][0]->copy());
		MITEM m01(m[0][1]->copy());
		MITEM m10(m[1][0]->copy());
		MITEM m11(m[1][1]->copy());
		MITEM D(matrix_determinant(m));
		if (D == 0.0) throw DivisionByZero();

		MMatrix mi;
		mi.Create(2,2);
		mi[0][0] = m11.copy();
		mi[1][1] = m00.copy();
		mi[0][1] = (-m01).copy();
		mi[1][0] = (-m10).copy();

		MMatrix& Mi = *(mi/D);
		return &Mi;
	}
	else if (nrows == 3)
	{
		MITEM m00(m[0][0]->copy());
		MITEM m01(m[0][1]->copy());
		MITEM m02(m[0][2]->copy());
		MITEM m10(m[1][0]->copy());
		MITEM m11(m[1][1]->copy());
		MITEM m12(m[1][2]->copy());
		MITEM m20(m[2][0]->copy());
		MITEM m21(m[2][1]->copy());
		MITEM m22(m[2][2]->copy());
		MITEM D(matrix_determinant(m));
		if (D == 0.0) throw DivisionByZero();

		MMatrix mi;
		mi.Create(3,3);
		mi[0][0] = (m11*m22 - m21*m12).copy();
		mi[1][0] = (m12*m20 - m10*m22).copy();
		mi[2][0] = (m10*m21 - m20*m11).copy();
		mi[0][1] = (m21*m02 - m01*m22).copy();
		mi[1][1] = (m00*m22 - m20*m02).copy();
		mi[2][1] = (m20*m01 - m00*m21).copy();
		mi[0][2] = (m01*m12 - m11*m02).copy();
		mi[1][2] = (m10*m02 - m00*m12).copy();
		mi[2][2] = (m00*m11 - m10*m01).copy();

		MMatrix& Mi  = *(mi/D);
		return &Mi;
	}
	else throw InvalidOperation();

	// we should not get here
	assert(false);
	return 0;
}
