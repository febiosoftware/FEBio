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
#include "math.h"
#include "fecore_api.h"

///////////////////////////////////////////////////////////////////////////////
// LINEAR SOLVER : colsol
// This solver solves linear system of matrices using a skyline format storage
// and a column reduction scheme. The symmetric matrix is stored in skyline format, 
// where the values array store the matrix elements that are below the skyline 
// and pointers is an array of indices that point to the diagonal elements of 
// the matrix.
// The matrix is overwritten with the LDLt factorization and the right hand side
// vector R is replaced by the solution.
//
// The implementation is split into two routines. A matrix factorization and a 
// back substitution part. colsol_factor performs the LDLt factorization while
// colsol_solve does the backsubstitution. colsol_solve needs to be called after
// colsol_factor. In order to solve for multiple right hand sides call colsol_factor
// once and then call colsol_solve with the different right hand sides.
//
// Details of the algorithm can be found in Bathe, "Finite Element Procedures",
// section 8.2, page 696 and following
//

FECORE_API void colsol_factor(int N, double* values, int* pointers)
{
	int i, j, r, mi, mj, mm;
	double krj;
	int pi, pj;

	// -A- factorize the matrix 

	// repeat over all columns
	for (j=1; j<N; ++j)
	{
		// find the first non-zero row in column j
		mj = j+1 - pointers[j+1] + pointers[j];

		pj = pointers[j]+j;

		// loop over all rows in column j
		for (i=mj+1; i<j; ++i)
		{
			// find the first non-zero row in column i
			mi = i+1 - pointers[i+1] + pointers[i];

			// determine max of mi and mj
			mm = (mi > mj ? mi : mj);

			pi = pointers[i]+i;

			double& kij = values[pj - i];

			// the next line is replaced by the piece of code between arrows
			// where the r loop is unrolled to give this algorithm a 
			// significant boost in speed. 
			// Although on good compilers this should not do much,
			// on compilers that do a poor optimization this trick can
			// double the speed of this algorithm.

//			for (r=mm; r<i; ++r) kij -= values[pi - r]*values[pj - r];

//-------------->
			for (r=mm; r<i-7; r+=8) 
			{
				kij -= values[pi - r  ]*values[pj - r  ] +
				       values[pi - r-1]*values[pj - r-1] +
				       values[pi - r-2]*values[pj - r-2] +
				       values[pi - r-3]*values[pj - r-3] +
				       values[pi - r-4]*values[pj - r-4] +
				       values[pi - r-5]*values[pj - r-5] +
				       values[pi - r-6]*values[pj - r-6] +
				       values[pi - r-7]*values[pj - r-7];
			}

			for (r=0; r<(i-mm)%8; ++r)
					kij -= values[pi - (i-1)+r]*values[pj - (i-1)+r];
//-------------->

		}

		// determine l[i][j]
		for (i=mj; i<j; ++i) values[pj - i] /= values[ pointers[i] ];

		// calculate d[j][j] value
		double& kjj = values[ pointers[j] ];
		for (r=mj; r<j; ++r) 
		{
			krj = values[pj - r];
			kjj -= krj*krj*values[ pointers[r] ];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

FECORE_API void colsol_solve(int N, double* values, int* pointers, double* R)
{
	int i, mi, r;

	// -B- back substitution

	// calculate V = L^(-T)*R vector
	for (i=1; i<N; ++i)
	{
		mi = i+1 - pointers[i+1] + pointers[i];
		for (r=mi; r<i; ++r)
			R[i] -= values[ pointers[i] + i - r]*R[r];
	}

	// calculate Vbar = D^(-1)*V
	for (i=0; i<N; ++i) R[i] /= values[ pointers[i] ];

	// calculate the solution
	for (i=N-1; i>0; --i)
	{
		mi = i+1 - pointers[i+1] + pointers[i];

		const double ri = R[i];
		const int pi = pointers[i] + i;

		// the following line was replaced
		// by the code segment between the arrows
//		for (r=mi; r<i; ++r) R[r] -= values[ pointers[i] + i - r]*R[i];

//--------->
		for (r=mi; r<i-7; r += 8) 
		{
			R[r  ] -= values[ pi - r  ]*ri;
			R[r+1] -= values[ pi - r-1]*ri;
			R[r+2] -= values[ pi - r-2]*ri;
			R[r+3] -= values[ pi - r-3]*ri;
			R[r+4] -= values[ pi - r-4]*ri;
			R[r+5] -= values[ pi - r-5]*ri;
			R[r+6] -= values[ pi - r-6]*ri;
			R[r+7] -= values[ pi - r-7]*ri;
		}

		for (r=0; r<(i-mi)%8; ++r)
			R[(i-1)- r] -= values[ pi - (i-1) + r  ]*ri;
//--------->
	}
}


///////////////////////////////////////////////////////////////////////////////
// This LU solver is grabbed from Numerical Recipes in C.
// To solve a system of equations first call ludcmp to calculate
// its LU decomposition. Next you can call the lubksb to perform
// a back substitution. The preferred way of using these functions
// therefore is:
//
// double** a;
// double* b;
// int* indx;
// int n;
// ...
// ludcmp(a,n,indx)
// lubksb(a,n,indx, b)
//
// The lubksb can be called as many times as there are rhs vectors for 
// this system that you wish to call. For a specific matrix ludcmp only
// has to be called once.
//

void ludcmp(double**a, int n, int* indx)
{
	int i, imax, j, k;
	double big, dum, sum, temp;
	double* vv;

	const double TINY = 1.0e-20;

	vv = new double[n];

	for (i=0; i<n; ++i)
	{
		big = 0;
		for (j=0; j<n; ++j)
			if ((temp = fabs(a[i][j])) > big) big = temp;

		vv[i] = 1.0 / big;
	}

	for (j=0; j<n; ++j)
	{
		for (i=0; i<j; ++i)
		{
			sum = a[i][j];
			for (k=0; k<i; ++k) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0;
		imax = j;
		for (i=j; i<n; ++i)
		{
			sum = a[i][j];
			for (k=0; k<j; ++k) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if ((dum=vv[i]*fabs(sum)) >= big)
			{
				big = dum;
				imax = i;
			}
		}
		if (j != imax)
		{
			for (k=0; k<n; ++k)
			{
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) 
		{
			a[j][j] = TINY;
		}
		if (j!=n-1)
		{
			dum = 1.0/(a[j][j]);
			for (i=j+1;i<n; ++i) a[i][j] *= dum;
		}
	}

	// clean up
	delete [] vv;
}

///////////////////////////////////////////////////////////////////////////////

void lubksb(double**a, int n, int *indx, double b[])
{
	int i, ii=0, ip, j;
	double sum;

	for (i=0; i<n; ++i)
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii != 0)
			for (j=ii-1; j<i; ++j) sum -= a[i][j]*b[j];
		else if (sum != 0.0) ii=i+1;
		b[i] = sum;
	}
	for (i=n-1; i>=0; --i)
	{
		sum = b[i];
		for (j=i+1;j<n;++j) sum -= a[i][j]*b[j];
		b[i] = sum/a[i][i];
	}
}
