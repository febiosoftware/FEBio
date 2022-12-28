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
#include "FEASTEigenSolver.h"
#include <FECore/CompactSymmMatrix.h>


BEGIN_FECORE_CLASS(FEASTEigenSolver, EigenSolver)
	ADD_PARAMETER(m_m0, "m0");
	ADD_PARAMETER(m_emin, "emin");
	ADD_PARAMETER(m_emax, "emax");
END_FECORE_CLASS();

#ifdef MKL_ISS
#undef PARDISO
#include <mkl.h>

FEASTEigenSolver::FEASTEigenSolver(FEModel* fem) : EigenSolver(fem)
{
	m_emin = 0.0;
	m_emax = 0.0;
	m_m0 = 0;
}

bool FEASTEigenSolver::Init()
{
	// initialize the FEAST solver
	feastinit(m_fpm);

	// fpm[0] = Specifies whether Extended Eigensolver routines print runtime status. (0 = no -default-, 1 = yes)
	// fpm[1] = Number of contour points Ne = 8. Must be one of {3,4,5,6,8,10,12,16,20,24,32,40,48}
	// fpm[2] = Error trace double precision stopping criteria eps = 10^-fpm[2] (default = 12)
	// fpm[3] = Maximum number of Extended Eigensolver refinement loops allowed. 
	//          If no convergence is reached within fpm[3] refinement loops, Extended Eigensolver routines return info=2
	// fpm[4] = User initial subspace. If fpm[4]=0 then Extended Eigensolver routines generate initial subspace, 
	//          if fpm[4]=1 the user supplied initial subspace is used.
	// fpm[5] = Extended Eigensolver stopping test.
	//          fpm[5]=0 (default): Extended Eigensolvers are stopped if this residual stopping test is satisfied.
	//          fpm[5]=1 : Extended Eigensolvers are stopped if this trace stopping test is satisfied.
	// fpm[6] = Error trace single precision stopping criteria(10 - fpm[6]) .(default = 5)
	// fpm[13] = fpm[13] = 0: Standard use for Extended Eigensolver routines. (default)
	//		     fpm[13] = 1: Non-standard use for Extended Eigensolver routines : return the computed eigenvectors subspace after one single contour integration.
	// fpm[26] = Specifies whether Extended Eigensolver routines check input matrices(applies to CSR format only). (0 = no -default-, 1 = yes)
	// fpm[27] = Check if matrix B is positive definite.Set fpm[27] = 1 to check if B is positive definite. (default = 0)
	// fpm[63] = Use the Intel MKL PARDISO solver with the user - defined PARDISO iparm array settings.
	//           fpm[63]=0 (default): Extended Eigensolver routines use the Intel MKL PARDISO default iparm settings defined by calling the pardisoinit subroutine.
	//	         fpm[63]=1 : The values from fpm[64] to fpm[127] correspond to iparm[0] to iparm[63] respectively according to the formula fpm[64 + i] = iparm[i] for i = 0, 1, ..., 63

#ifdef _DEBUG
	m_fpm[0] = 1; // turn on FEAST output
#endif

	return true;
}

bool FEASTEigenSolver::EigenSolve(SparseMatrix* A, SparseMatrix* B, vector<double>& eigenValues, matrix& eigenVectors)
{
	CompactSymmMatrix* cmA = dynamic_cast<CompactSymmMatrix*>(A);
	if (cmA == nullptr) return false;

	CompactSymmMatrix* cmB = nullptr; 
	if (B)
	{
		cmB = dynamic_cast<CompactSymmMatrix*>(B);
		if (cmB == nullptr) return false;

		// make sure the matrices have the same size
		if (cmA->Rows() != cmB->Rows()) return false;
	}

	// parameters to dfeast_scsrgv
	// Input parameters:
	const char* uplo = "U"; // the matrix is stored as lower triangular. 
	MKL_INT n = cmA->Rows();
	double* a = cmA->Values(); // nonzero values of matrix A
	MKL_INT* ai = cmA->Pointers();
	MKL_INT* aj = cmA->Indices(); 

	// Output parameters
	MKL_INT m = 0; // number of eigenvalues found in [emin, emax], 0 <= m <= m0.
	double epsout = 0; // relative error on the trace
	MKL_INT loops = 0; // number of refinement loops executed.

	eigenValues.resize(m_m0, 0.0);
	double* e = &eigenValues[0];

	eigenVectors.resize(n, m_m0);
	double* x = eigenVectors[0];

	vector<double> res(m_m0); // the first m entries contain the relative residual error. 

	MKL_INT info = 0; // return value. If 0, all is well. 

	if (B)
	{
		double* b = cmB->Values(); // nonzero values of matrix B
		MKL_INT* bi = cmB->Pointers();
		MKL_INT* bj = cmB->Indices();
		dfeast_scsrgv(uplo, &n, a, ai, aj, b, bi, bj, m_fpm, &epsout, &loops, &m_emin, &m_emax, &m_m0, e, x, &m, &res[0], &info);
	}
	else
	{
		dfeast_scsrev(uplo, &n, a, ai, aj, m_fpm, &epsout, &loops, &m_emin, &m_emax, &m_m0, e, x, &m, &res[0], &info);
	}

	if (m == 0)
	{
		eigenValues.clear();
	}
	else if (m != m_m0) eigenValues.resize(m);

	return (info == 0);
}
#else 
FEASTEigenSolver::FEASTEigenSolver(FEModel* fem) : EigenSolver(fem)  {}
bool FEASTEigenSolver::Init() { return false; }
bool FEASTEigenSolver::EigenSolve(SparseMatrix* A, SparseMatrix* B, vector<double>& eigenValues, matrix& eigenVectors) { return false; }
#endif