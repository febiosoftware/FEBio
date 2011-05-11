#include "stdafx.h"
#include <math.h>

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
