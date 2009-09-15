#include "stdafx.h"
#include "LUSolver.h"

//-----------------------------------------------------------------------------
bool LUSolver::PreProcess(SparseMatrix& K)
{
	FullMatrix* pK = dynamic_cast<FullMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix needs to be a FullMatrix for this solver\n");
		return false;
	}

	// We don't need to do any preprocessing for this solver

	return LinearSolver::PreProcess(K);
}

//-----------------------------------------------------------------------------
bool LUSolver::Factor(SparseMatrix& K)
{
	// convert to a FullMatrix
	FullMatrix& a = dynamic_cast<FullMatrix&> (K);

	const double TINY = 1.0e-20;
	int i, imax, j, k;
	double big, dum, sum, temp;

	int n = a.Size();
	// create index vector
	indx.create(n);

	vector<double> vv(n);
	for (i=0; i<n; ++i)
	{
		big = 0;
		for (j=0; j<n; ++j)
			if ((temp=fabs(a(i,j))) > big) big = temp;
		if (big == 0) return false; // singular matrix
		vv[i] = 1.0 / big;
	}

	for (j=0; j<n; ++j)
	{
		for (i=0; i<j; ++i)
		{
			sum = a(i,j);
			for (k=0; k<i; ++k) sum -= a(i,k)*a(k,j);
			a(i,j) = sum;
		}
		big = 0;
		imax = j;
		for (i=j;i<n;++i)
		{
			sum = a(i,j);
			for (k=0; k<j; ++k) sum -= a(i,k)*a(k,j);
			a(i,j) = sum;
			if ((dum=vv[i]*fabs(sum))>=big)
			{
				big = dum;
				imax = i;
			}
		}

		if (j != imax)
		{
			for (k=0; k<n; ++k)
			{
				dum = a(imax,k);
				a(imax,k) = a(j,k);
				a(j,k) = dum;
			}
			vv[imax] = vv[j];
		}

		indx[j] = imax;
		if (a(j,j) == 0) a(j,j) = TINY;
		if (j != n-1)
		{
			dum = 1.0/a(j,j);
			for (i=j+1;i<n; ++i) a(i,j) *= dum;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool LUSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& b)
{
	FullMatrix& a = dynamic_cast<FullMatrix&> (K);

	x = b;

	int i, ii=0, ip, j;
	double sum;

	int n = a.Size();
	for (i=0; i<n; ++i)
	{
		ip = indx[i];
		sum = x[ip];
		x[ip] = x[i];
		if (ii != 0)
			for (j=ii-1;j<i;++j) sum -= a(i,j)*x[j];
		else if (sum != 0)
			ii = i+1;
		x[i] = sum;
	}

	for (i=n-1; i>=0; --i)
	{
		sum = x[i];
		for (j=i+1; j<n; ++j) sum -= a(i,j)*x[j];
		x[i] = sum/a(i,i);
	}

	return false;
}

//-----------------------------------------------------------------------------
bool LUSolver::Solve(SparseMatrix& K, matrix& x, matrix& b)
{
	return false;
}

//-----------------------------------------------------------------------------
void LUSolver::Destroy(SparseMatrix& K)
{
	// nothing to destroy
	LinearSolver::Destroy(K);
}
