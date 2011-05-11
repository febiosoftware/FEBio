#include "stdafx.h"


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

void colsol_factor(int N, double* values, int* pointers)
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

void colsol_solve(int N, double* values, int* pointers, double* R)
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
