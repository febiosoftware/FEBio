// LinSolver.cpp: implementation of the LinSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include "LinearSolver.h"

void colsol_factor(int N, double* values, int* pointers);
void colsol_solve (int N, double* values, int* pointers, double* R);

//////////////////////////////////////////////////////////////////////
// LUSolver
//////////////////////////////////////////////////////////////////////

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

bool LUSolver::Solve(SparseMatrix& K, matrix& x, matrix& b)
{
	return false;
}

void LUSolver::Destroy()
{
	// nothing to destroy
	LinearSolver::Destroy();
}

//////////////////////////////////////////////////////////////////////
// SkylineSolver
//////////////////////////////////////////////////////////////////////

bool SkylineSolver::PreProcess(SparseMatrix& K)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n");
		return false;
	}

	// We don't need to do any preprocessing for this solver

	return LinearSolver::PreProcess(K);
}

bool SkylineSolver::Factor(SparseMatrix& K)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n");
		return false;
	}
#ifdef DEBUG
	int i, n;
	int* pointers;
	double* values;

	n = pK->Size();
	pointers = pK->pointers();
	values = pK->values();

	printf("\nPointers, Values");
	for (i=0; i<n; i++) printf("\n%d, %g", pointers[i], values[i]);
#endif

	colsol_factor(pK->Size(), pK->values(), pK->pointers());

	return true;
}

bool SkylineSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& R)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n");
		return false;
	}

	// we need to make a copy of R since colsol overwrites the right hand side vector
	// with the solution
	int neq = pK->Size();
	for (int i=0; i<neq; ++i) x[i] = R[i];
	colsol_solve(pK->Size(), pK->values(), pK->pointers(), x);

	return true;
}

bool SkylineSolver::Solve(SparseMatrix& K, matrix& x, matrix& R)
{
	// Let's make sure the matrix K is of the correct type
	SkylineMatrix* pK = dynamic_cast<SkylineMatrix*> (&K);

	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n");
		return false;
	}

	// we need to make a copy of R since colsol overwrites the right hand side vector
	// with the solution
	int nrhs = x.rows();
	for (int i=0; i<nrhs; ++i)
	{
		int neq = pK->Size();
		for (int j=0; j<neq; ++j) x[i][j] = R[i][j];

		colsol_solve(K.Size(), pK->values(), pK->pointers(), x[i]);
	}

	return true;
}


void SkylineSolver::Destroy()
{
	// Nothing to destroy
	LinearSolver::Destroy();
}

//////////////////////////////////////////////////////////////////////
// PSLDLTSolver
//////////////////////////////////////////////////////////////////////

bool PSLDLTSolver::PreProcess(SparseMatrix& K)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// let's make sure the matrix K is of the correct type
	CompactMatrix* pK = dynamic_cast<CompactMatrix*> (&K);
	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n\n");
		return false;
	}

	// Do the preprocessing
	int nonz;
	double ops;
	PSLDLT_Preprocess(0, pK->Size(), pK->pointers(), pK->indices(), &nonz, &ops);

	return LinearSolver::PreProcess(K);

#endif

}

bool PSLDLTSolver::Factor(SparseMatrix& K)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// let's make sure the matrix K is of the correct type
	CompactMatrix* pK = dynamic_cast<CompactMatrix*> (&K);
	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n\n");
		return false;
	}

#ifdef DEBUG
	printf("\nThis is a test");

	int i, n, nnz, *pointers, *indices;
	double* values;

	fflush(stdout);
	n = pK->Size();
	nnz = pK->NonZeroes();
	pointers = pK->pointers();
	indices = pK->indices();
	values = pK->values();
	fprintf(stdout, "\nPointers:");
	for (i=0; i<n; i++) fprintf(stdout, "\n%d", pointers[i]);
	fprintf(stdout, "\nIndices, Values:");
	for (i=0; i<nnz; i++) fprintf(stdout, "\n%d, %g", indices[i], values[i]);
#endif

	// Do the factorization
	PSLDLT_Factor(0, pK->Size(), pK->pointers(), pK->indices(), pK->values());
	return true;

#endif

}

bool PSLDLTSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& R)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// let's make sure the matrix K is of the correct type
	CompactMatrix* pK = dynamic_cast<CompactMatrix*> (&K);
	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n\n");
		return false;
	}

	// Let's roll !!
	PSLDLT_Solve(0, x, R);

	return true;

#endif
}

bool PSLDLTSolver::Solve(SparseMatrix& K, matrix& x, matrix& b)
{
	// First, make sure the PSLDLT solver is available on this platform
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
	return false;
#else

	// let's make sure the matrix K is of the correct type
	CompactMatrix* pK = dynamic_cast<CompactMatrix*> (&K);
	if (pK == 0)
	{
		fprintf(stderr, "Stiffness matrix is not of correct type for this solver\n\n");
		return false;
	}


	// Let's roll !!
	int nrhs = x.rows();
	for (int i=0; i<nrhs; ++i) PSLDLT_Solve(0, x[i], b[i]);

	return true;
#endif
}

void PSLDLTSolver::Destroy()
{
#ifndef PSLDLT
	fprintf(stderr, "FATAL ERROR : The PSLDLT solver is not available on this platform\n\n");
#else
	if (m_bvalid) PSLDLT_Destroy(0);
	LinearSolver::Destroy();
#endif
}

//////////////////////////////////////////////////////////////////////
// WSMPSolver
//////////////////////////////////////////////////////////////////////

bool WSMPSolver::PreProcess(SparseMatrix& K)
{
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform.\n\n");
	return false;
#else

	return false;
#endif
}

bool WSMPSolver::Factor(SparseMatrix& K)
{
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform.\n\n");
	return false;
#else

	return false;
#endif
}

bool WSMPSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& b)
{
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform.\n\n");
	return false;
#else

	return false;
#endif
}

bool WSMPSolver::Solve(SparseMatrix& K, matrix& x, matrix& b)
{
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform.\n\n");
	return false;
#else

	return false;
#endif
}

void WSMPSolver::Destroy()
{
#ifndef WSMP
	fprintf(stderr, "FATAL ERROR: The WSMP solver is not available on this platform.\n\n");
	return;
#else

#endif
}

//-----------------------------------------------------------------------------

ConjGradIterSolver::ConjGradIterSolver()
{
	m_tol = 0.01;
	m_kmax = 200;
	m_nprint = 0;
}

bool ConjGradIterSolver::PreProcess(SparseMatrix& K)
{
	return true;
}

bool ConjGradIterSolver::Factor(SparseMatrix& K)
{
/*	CompactMatrix& C = dynamic_cast<CompactMatrix&>(K);
	int neq = K.Size();

	int* pi = C.indices();
	int* pr = C.pointers();
	double* pv = C.values();

	// copy the diagonals to P
	m_P.create(neq);
	int i, j, n, l;

	// store the diagonal elements in P
	for (j=0; j<neq; ++j) m_P[j] = 1.0 / pv[pr[j]];

	for (j=0; j<neq; ++j)
	{
		// get the number of entries on the column
		l = pr[j+1]-pr[j];

		// loop over all rows
		for (i=0; i<l; ++i)
		{
			n = pi[pr[j] + i];
			pv[ pr[j] + i] *= m_P[n];
		}
	}
*/
	return true;
}

bool ConjGradIterSolver::Solve(SparseMatrix& K, matrix& x, matrix& b)
{
	assert(false);
	return false;
}

bool ConjGradIterSolver::Solve(SparseMatrix& K, vector<double>& x, vector<double>& b)
{
	int i;

	CompactMatrix& A = dynamic_cast<CompactMatrix&>(K);

	int N = x.size();

	// intialize x to zero
	x.zero();

	// iteration counter
	int k = 0;

	// initial residual vector
	vector<double> r(N);
	A.mult_vector(x, r);
	for (i=0; i<N; ++i) r[i] = b[i] - r[i];

	// initial norms
	double rho1 = r*r, rho2;
	double normb = 0;
	for (i=0; i<N; ++i) normb += (b[i])*(b[i]);
	normb = sqrt(normb);

	// search direction
	vector<double> p(N);

	// work vector
	vector<double> w(N);

	double alpha, beta;

	if (m_nprint > 0) fprintf(stderr, "Solving linear system ...\n");

	// loop until converged or until max iterations reached
	while ((sqrt(rho1) > m_tol*normb) && (k < m_kmax))
	{
		// increase iteration counter
		++k;

		// calculate search direction
		if (k==1)
		{
			p = r;
		}
		else
		{
			beta = rho1/rho2;
			for (i=0; i<N; ++i) p[i] = r[i] + p[i]*beta;
		}
		A.mult_vector(p, w);
		alpha = rho1 / (p*w);

		for (i=0; i<N; ++i)
		{
			x[i] += p[i]*alpha;
			r[i] -= w[i]*alpha;
		}
		rho2 = rho1;
		rho1 = r*r;

		if (m_nprint > 0)
			fprintf(stderr, "%d: %lg\n", k, sqrt(rho1)/normb);
	}

	if (m_nprint > 0)
	{
		if ((k >= m_kmax) && (sqrt(rho1) > m_tol*normb))
		{
			fprintf(stderr, "Max iterations reached. Solution has not converged.\n");
		}
		else fprintf(stderr, "Solution converged\n");
	}

	return true;
}

void ConjGradIterSolver::Destroy()
{

}
