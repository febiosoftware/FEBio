// MatrixFactory.cpp: implementation of the MatrixFactory class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MatrixFactory.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MatrixFactory::MatrixFactory()
{

}

MatrixFactory::~MatrixFactory()
{

}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : MatrixFactory::CreateMatrix
// This function allocates the initial storage for the spare matrix A. RTTI is 
// used the resolve the exact type of A. Based in its type the corresponding
// creation routine is invoked.
// The LM parameter contains for each element the global equation numbers for the
// element's degrees of freedom. 'neq' is the number of equations. 
//

bool MatrixFactory::CreateMatrix(SparseMatrix* pA, vector< vector<int> >& LM, int neq)
{
	// find out what type of matrix we're dealing with
	CompactUnSymmMatrix* pCompactUS;
	CompactMatrix* pCompact;
	SkylineMatrix* pSkyline;
	FullMatrix* pFull;

	if      (pSkyline   = dynamic_cast<SkylineMatrix*      > (pA)) return CreateSkyline      (pSkyline  , LM, neq);
	else if (pCompact   = dynamic_cast<CompactMatrix*      > (pA)) return CreateCompact      (pCompact  , LM, neq);
	else if (pCompactUS = dynamic_cast<CompactUnSymmMatrix*> (pA)) return CreateCompactUnSymm(pCompactUS, LM, neq);
	else if (pFull      = dynamic_cast<FullMatrix*         > (pA)) return CreateFull         (pFull     , LM, neq);
	
	// If we get here that means the matrix pA is not of a known type.
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : MatrixFactory::CreateMatrix
// This function allocates the initial storage for the spare matrix A. RTTI is 
// used the resolve the exact type of A. Based in its type the corresponding
// creation routine is invoked.
// The LM parameter contains for each element the global equation numbers for the
// element's degrees of freedom. 'neq' is the number of equations. 
//

bool MatrixFactory::CreateMatrix(SparseMatrix* pA, MatrixProfile& mp)
{
	// find out what type of matrix we're dealing with
	CompactUnSymmMatrix* pCompactUS;
	CompactMatrix* pCompact;
	SkylineMatrix* pSkyline;
	FullMatrix* pFull;

	if      (pSkyline   = dynamic_cast<SkylineMatrix*      > (pA)) return CreateSkyline      (pSkyline  , mp);
	else if (pCompact   = dynamic_cast<CompactMatrix*      > (pA)) return CreateCompact      (pCompact  , mp);
	else if (pCompactUS = dynamic_cast<CompactUnSymmMatrix*> (pA)) return CreateCompactUnSymm(pCompactUS, mp);
	else if (pFull      = dynamic_cast<FullMatrix*         > (pA)) return CreateFull         (pFull     , mp);
	
	// If we get here that means the matrix pA is not of a known type.
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : MatrixFactory::CreateFull
// This routine allocates storage for a full matrix.
//

bool MatrixFactory::CreateFull(FullMatrix* pA, vector< vector<int> >& LM, int neq)
{
	pA->Create(neq);
	return true;
}

bool MatrixFactory::CreateFull(FullMatrix* pA, MatrixProfile& mp)
{
	pA->Create(mp.size());
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : MatrixFactory::CreateSkyline
// This routine allocates storage for a skyline matrix.
//

bool MatrixFactory::CreateSkyline(SkylineMatrix* pA, vector< vector<int> >& LM, int neq)
{
	int i, j, k, N;

	int* pointers = new int[neq+1];
	if (pointers == 0) throw MemException(sizeof(int)*(neq+1));

	for (i=0; i<neq; ++i) pointers[i] = i;

	int* lm;

	int M = LM.size();

	// repeat for all elements
	for (i=0; i<M; ++i)
	{
		lm = LM[i];
		N = LM[i].size();
	
		for (j=0; j<N; ++j)
			if (lm[j] >= 0)
			{
				for (k=0; k<N; ++k)
				{
					if (lm[k] >= 0 && lm[k] < pointers[lm[j]])
						pointers[ lm[j] ] = lm[k];
				}
			}
	}

	// setup pointers array
	for (i=0; i<neq; ++i) pointers[i] = i - pointers[i]+1;
	for (i=1; i<neq; ++i) pointers[i] += pointers[i-1];
	for (i=neq; i>0; --i) pointers[i] = pointers[i-1];
	pointers[0] = 0;

	// allocate stiffness matrix
	double* values = new double[pointers[neq]];
	if (values==0) throw MemException(sizeof(double)*neq);

	// create the matrix
	pA->Create(values, pointers, neq);

	return true;
}

bool MatrixFactory::CreateSkyline(SkylineMatrix* pA, MatrixProfile& mp)
{
	int i, n;

	int neq = mp.size();
	int* pointers = new int[neq + 1];
	if (pointers == 0) throw MemException(sizeof(int)*(neq+1));

	pointers[0] = 0;
	for (i=1; i<=neq; ++i)
	{
		vector<int>& a = mp.column(i-1);
		n = i - a[0];
		pointers[i] = pointers[i-1] + n;
	}

	// allocate stiffness matrix
	double* values = new double[pointers[neq]];
	if (values==0) 
	{
		double falloc = (double) sizeof(double) * (double) (pointers[neq]);
		throw MemException(falloc);
	}

	// create the matrix
	pA->Create(values, pointers, neq);

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : MatrixFactory::CreateCompact
// This routine allocates storage for a compact matrix (Hartwell-Boeing format)
//

bool MatrixFactory::CreateCompact(CompactMatrix* pA, vector< vector<int> >& LM, int neq)
{
	int i, j, k, iel, n, *lm, N, Ntot = 0;

	int M = LM.size();

	// Find the number of elements that contribute to a certain column
	int* pval = new int[neq];
	if (pval == 0) throw MemException(sizeof(int)*neq);

	for (i=0; i<neq; ++i) pval[i] = 0;

	for (i=0; i<M; ++i)
	{
		lm = LM[i];
		N = LM[i].size();
		Ntot += N;
		for (j=0; j<N; ++j)
		{
			if (lm[j] >= 0) pval[lm[j]]++; 
		}
	}

	// create a "compact" 2D array that stores for each column the element
	// numbers that contribute to that column
	int* pelc = new int[Ntot];
	if (pelc == 0) throw MemException(sizeof(int)*Ntot);

	int** ppelc = new int*[neq];
	if (ppelc == 0) throw MemException(sizeof(int*)*neq);

	ppelc[0] = pelc;
	for (i=1; i<neq; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	for (i=0; i<M; ++i)
	{
		lm = LM[i];
		N = LM[i].size();
		for (j=0; j<N; ++j)
		{
			if (lm[j] >= 0) *(ppelc[lm[j]])++ = i;
		}
	}

	// reset pelc pointers
	ppelc[0] = pelc;
	for (i=1; i<neq; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	// create the pointers array
	// we use the pointers array to store the size of the columns
	// note that we store these values shifted up one column
	// we do this so it will be easier to correct the pointers array later
	int* pointers = new int[neq + 1];
	if (pointers == 0) throw MemException(sizeof(int)*(neq+1));

	pointers[0] = 0;	

	// loop over all columns, flaggin which row items are used
	// keeping track of the total count as well
	int* pcol = new int[neq];
	if (pcol == 0) throw MemException(sizeof(int)*neq);

	int nsize = 0; // --> size of indices and values array
	for (i=0; i<neq; ++i)
	{
		// zero column
		for (j=i; j<neq; ++j) pcol[j] = 0;
		n = 0;
	
		// loop over all elements in the plec
		for (j=0; j<pval[i]; ++j)
		{
			iel = (ppelc[i])[j];
			lm = LM[iel];
			N = LM[iel].size();
			for (k=0; k<N; ++k)
			{
				if ((lm[k] >= i) && (pcol[ lm[k] ] == 0)) { n++; pcol[ lm[k] ] = 1; }
			}
		}

		pointers[i+1] = pointers[i] + n;
		nsize += n;
	}

	// create the indices array
	int* pindices = new int[nsize];
	if (pindices == 0) throw MemException(sizeof(int)*nsize);

	int* pi;

	for (i=0; i<neq; ++i)
	{
		// zero column
		for (j=i; j<neq; ++j) pcol[j] = 0;

		// loop over all elements in the plec
		for (j=0; j<pval[i]; ++j)
		{
			iel = ppelc[i][j];
			lm = LM[iel];
			N = LM[iel].size();
			for (k=0; k<N; ++k)
			{
				if ((lm[k] >= i) && (pcol[ lm[k] ] == 0)) { pcol[ lm[k] ] = 1; }
			}
		}

		// we do the following step to obtain a sorted column list
		n = 0;
		pi = pindices + pointers[i];
		for (j=i; j<neq; ++j)
		{
			if (pcol[j] == 1) { pi[n] = j; ++n; }
		}
	}

	// cleanup
	delete [] pcol;
	delete [] ppelc;
	delete [] pelc;
	delete [] pval;

	// create the values array
	double* pvalues = new double[nsize];
	if (pvalues == 0) throw MemException(sizeof(double)*nsize);

	// create the stiffness matrix
	pA->Create(neq, nsize, pvalues, pindices, pointers);

	return true;
}

bool MatrixFactory::CreateCompact(CompactMatrix* pA, MatrixProfile& mp)
{
	int i, j, k, n;

	int neq = mp.size();

	int* pointers = new int[neq + 1];
	for (i=0; i<=neq; ++i) pointers[i] = 0;

	int nsize = 0;
	for (i=0; i<neq; ++i)
	{
		vector<int>& a = mp.column(i);
		n = a.size();
		for (j=0; j<n; j += 2)
		{
			nsize += a[j+1] - a[j] + 1;
			for (k=a[j]; k<=a[j+1]; ++k) pointers[k]++;
		}
	}

	int* pindices = new int[nsize];
	int m = 0;
	for (i=0; i<=neq; ++i)
	{
		n = pointers[i];
		pointers[i] = m;
		m += n;
	}

	int* pval = new int[neq];
	for (i=0; i<neq; ++i) pval[i] = 0;

	for (i=0; i<neq; ++i)
	{
		vector<int>& a = mp.column(i);
		n = a.size();
		for (j=0; j<n; j += 2)
		{
			for (k=a[j]; k<=a[j+1]; ++k) 
			{
				pindices[ pointers[k] + pval[k]] = i;
				++pval[k];
			}
		}
	}

	// cleanup
	delete [] pval;

	// create the values array
	double* pvalues = new double[nsize];
	if (pvalues == 0) throw MemException(sizeof(double)*nsize);

	// create the stiffness matrix
	pA->Create(neq, nsize, pvalues, pindices, pointers);


	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : MatrixFactory::CreateCompact
// This routine allocates storage for a compact matrix (Hartwell-Boeing format)
//

bool MatrixFactory::CreateCompactUnSymm(CompactUnSymmMatrix* pA, vector< vector<int> >& LM, int neq)
{
	int i, j, k, iel, n, *pn, N, Ntot = 0;

	int M = LM.size();

	// Find the number of elements that contribute to a certain column
	int* pval = new int[neq];
	if (pval==0) throw MemException(sizeof(int)*neq);

	for (i=0; i<neq; ++i) pval[i] = 0;

	for (i=0; i<M; ++i)
	{
		pn = LM[i];
		N = LM[i].size();
		Ntot += N;
		for (j=0; j<N; ++j)
		{
			if (pn[j] >= 0) pval[pn[j]]++; 
		}
	} 

	// create a "compact" 2D array that stores for each column the element
	// numbers that contribute to that column
	int* pelc = new int[Ntot];
	if (pelc==0) throw MemException(sizeof(int)*Ntot);

	int** ppelc = new int*[neq];
	if (ppelc==0) throw MemException(sizeof(int*)*neq);

	ppelc[0] = pelc;
	for (i=1; i<neq; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	for (i=0; i<M; ++i)
	{
		pn = LM[i];
		N = LM[i].size();
		for (j=0; j<N; ++j)
		{
			if (pn[j] >= 0) *(ppelc[pn[j]])++ = i;
		}
	}

	// reset pelc pointers
	ppelc[0] = pelc;
	for (i=1; i<neq; ++i) ppelc[i] = ppelc[i-1] + pval[i-1];

	// create the pointers array
	// we use the pointers array to store the size of the columns
	// note that we store these values shifted up one column
	// we do this so it will be easier to correct the pointers array later
	int* pointers = new int[neq + 1];
	if (pointers==0) throw MemException(sizeof(int)*(neq+1));

	pointers[0] = 0;	

	// loop over all columns, flaggin which row items are used
	// keeping track of the total count as well
	int* pcol = new int[neq];
	if (pcol==0) throw MemException(sizeof(int)*neq);

	int nsize = 0; // --> size of indices and values array
	for (i=0; i<neq; ++i)
	{
		// zero column
		for (j=0; j<neq; ++j) pcol[j] = 0;
		n = 0;
	
		// loop over all elements in the plec
		for (j=0; j<pval[i]; ++j)
		{
			iel = (ppelc[i])[j];
			pn = LM[iel];
			N = LM[iel].size();
			for (k=0; k<N; ++k)
			{
				if (pcol[ pn[k] ] == 0) { n++; pcol[ pn[k] ] = 1; }
			}
		}

		pointers[i+1] = n;
		nsize += n;
	}

	// adjust pointers array
	for (i=0; i<neq; ++i) pointers[i+1] += pointers[i];

	// create the indices array
	int* pindices = new int[nsize];
	if (pindices==0) throw MemException(sizeof(int)*nsize);

	int* pi;

	for (i=0; i<neq; ++i)
	{
		// zero column
		for (j=0; j<neq; ++j) pcol[j] = 0;

		n = 0;
		pi = pindices + pointers[i];
	
		// loop over all elements in the plec
		for (j=0; j<pval[i]; ++j)
		{
			iel = ppelc[i][j];
			pn = LM[iel];
			N = LM[iel].size();
			for (k=0; k<N; ++k)
			{
				if (pcol[ pn[k] ] == 0) { pi[n] = pn[k]; n++; pcol[ pn[k] ] = 1; }
			}
		}
	}

	// cleanup
	delete [] pval;
	delete [] pelc;
	delete [] ppelc;
	delete [] pcol;

	// create the values array
	double* pvalues = new double[nsize];
	if (pvalues==0) throw MemException(sizeof(double)*nsize);

	// create the stiffness matrix
	pA->Create(neq, nsize, pvalues, pindices, pointers);

	return true;
}

bool MatrixFactory::CreateCompactUnSymm(CompactUnSymmMatrix* pA, MatrixProfile& mp)
{
	int i, j, k, n;

	int neq = mp.size();

	int* pointers = new int[neq + 1];
	for (i=0; i<=neq; ++i) pointers[i] = 0;

	int nsize = 0;
	for (i=0; i<neq; ++i)
	{
		vector<int>& a = mp.column(i);
		n = a.size();
		for (j=0; j<n; j += 2)
		{
			nsize += 2*(a[j+1] - a[j] + 1);
			pointers[i] += a[j+1] - a[j] + 1;
			for (k=a[j]; k<=a[j+1]; ++k) pointers[k]++;
		}
		--pointers[i]; // we double counted the diagonal
		--nsize;
	}

	int* pindices = new int[nsize];
	int m = 0;
	for (i=0; i<=neq; ++i)
	{
		n = pointers[i];
		pointers[i] = m;
		m += n;
	}
	assert(pointers[neq] == nsize);

	int* pval = new int[neq];
	for (i=0; i<neq; ++i) pval[i] = 0;

	for (i=0; i<neq; ++i)
	{
		vector<int>& a = mp.column(i);
		n = a.size();
		for (j=0; j<n; j += 2)
		{
			for (k=a[j]; k<=a[j+1]; ++k) 
			{
				pindices[ pointers[i] + pval[i]] = k;
				++pval[i];
			}
		}
		for (j=0; j<n; j += 2)
		{
			for (k=a[j]; k<=a[j+1]; ++k) 
			{
				if (k != i)
				{
					pindices[ pointers[k] + pval[k]] = i;
					++pval[k];
				}
			}
		}
	}

	// cleanup
	delete [] pval;

	// create the values array
	double* pvalues = new double[nsize];
	if (pvalues == 0) throw MemException(sizeof(double)*nsize);

	// create the stiffness matrix
	pA->Create(neq, nsize, pvalues, pindices, pointers);

	return true;
}

//-----------------------------------------------------------------------------
// FUNCTION : MatrixFactory::Assemble
// This function determines the explicit type of sparse matrix and calls
// the correct assembly routine
//
void MatrixFactory::Assemble(SparseMatrix& K, matrix& ke, vector<int>& LM)
{
	// find out what type of matrix we're dealing with
	CompactUnSymmMatrix* pCompactUS;
	CompactMatrix* pCompact;
	SkylineMatrix* pSkyline;
	FullMatrix* pFull;

	if      (pCompact   = dynamic_cast<CompactMatrix*      > (&K)) AssembleCompact      (*pCompact  , ke, LM);
	else if (pCompactUS = dynamic_cast<CompactUnSymmMatrix*> (&K)) AssembleCompactUnSymm(*pCompactUS, ke, LM);
	else if (pSkyline   = dynamic_cast<SkylineMatrix*      > (&K)) AssembleSkyline      (*pSkyline  , ke, LM);
	else if (pFull      = dynamic_cast<FullMatrix*         > (&K)) AssembleFull         (*pFull     , ke, LM);
}

//-----------------------------------------------------------------------------
// FUNCTION : MatrixFactory::Assemble
// This function determines the explicit type of sparse matrix and calls
// the correct assembly routine
//
void MatrixFactory::Assemble(SparseMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	// find out what type of matrix we're dealing with
	CompactUnSymmMatrix* pCompactUS;
	CompactMatrix* pCompact;
	SkylineMatrix* pSkyline;
	FullMatrix* pFull;

	if      (pCompact   = dynamic_cast<CompactMatrix*      > (&K)) AssembleCompact      (*pCompact  , ke, LMi, LMj);
	else if (pCompactUS = dynamic_cast<CompactUnSymmMatrix*> (&K)) AssembleCompactUnSymm(*pCompactUS, ke, LMi, LMj);
	else if (pSkyline   = dynamic_cast<SkylineMatrix*      > (&K)) AssembleSkyline      (*pSkyline  , ke, LMi, LMj);
	else if (pFull      = dynamic_cast<FullMatrix*         > (&K)) AssembleFull         (*pFull     , ke, LMi, LMj);
}

//-----------------------------------------------------------------------------
// FUNCTION : MatrixFactory::AssembleFull
// This function assembles the local stiffness matrix
// into the global stiffness matrix which is in full format
//

void MatrixFactory::AssembleFull(FullMatrix& K, matrix& ke, vector<int>& LM)
{
	int i, j, I, J;

	const int N = ke.rows();

	for (i=0; i<N; ++i)
	{
		if ((I = LM[i])>=0)
		{
			for (j=0; j<N; ++j)
			{
				if ((J = LM[j]) >= 0) K(I,J) += ke[i][j];
			}
		}
	}
}

void MatrixFactory::AssembleFull(FullMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	int i, j, I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	for (i=0; i<N; ++i)
	{
		if ((I = LMi[i])>=0)
		{
			for (j=0; j<M; ++j)
			{
				if ((J = LMj[j]) >= 0) K(I,J) += ke[i][j];
			}
		}
	}

}


//-----------------------------------------------------------------------------
// FUNCTION : MatrixFactory::AssembleSkyline
// This function assembles the local stiffness matrix
// into the global stiffness matrix which is in skyline format
//
void MatrixFactory::AssembleSkyline(SkylineMatrix& K, matrix& ke, vector<int>& LM)
{
	int i, j, I, J;

	const int N = ke.rows();

	double* pv = K.values();
	int* pi = K.pointers();

	for (i=0; i<N; ++i)
	{
		I = LM[i];

		if (I>=0)
		{
			for (j=0; j<N; ++j)
			{
				J = LM[j];

				// only add values to upper-diagonal part of stiffness matrix
				if (J>=I)
				{
					pv[ pi[J] + J - I] += ke[i][j];
				}
			}
		}
	}
}

void MatrixFactory::AssembleSkyline(SkylineMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	int i, j, I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	double* pv = K.values();
	int* pi = K.pointers();

	for (i=0; i<N; ++i)
	{
		I = LMi[i];

		if (I>=0)
		{
			for (j=0; j<M; ++j)
			{
				J = LMj[j];

				// only add values to upper-diagonal part of stiffness matrix
				if (J>=I)
				{
					pv[ pi[J] + J - I] += ke[i][j];
				}
			}
		}
	}

}

//-----------------------------------------------------------------------------
// FUNCTION : MatrixFactory::AssembleCompact
// This function assembles the local stiffness matrix
// into the global stiffness matrix which is in compact column storage
//
void MatrixFactory::AssembleCompact(CompactMatrix& K, matrix& ke, vector<int>& LM)
{
	int i, j, I, J;

	const int N = ke.rows();

	int* indices = K.indices();
	int* pointers = K.pointers();
	double* pd = K.values();

	int *pi, l, n;

	for (i=0; i<N; ++i)
	{
		I = LM[i];

		for (j=0; j<N; ++j)
		{
			J = LM[j];

			// only add values to lower-diagonal part of stiffness matrix
			if ((I>=J) && (J>=0)) 
			{
				pi = indices + pointers[J];
				l = pointers[J+1] - pointers[J];
				for (n=0; n<l; ++n) if (pi[n] == I) 
				{
					pd[pointers[J] + n] += ke[i][j];
					break;
				}
			}
		}
	}
}

void MatrixFactory::AssembleCompact(CompactMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	int i, j, I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	int* indices = K.indices();
	int* pointers = K.pointers();
	double* pd = K.values();

	int *pi, l, n;

	for (i=0; i<N; ++i)
	{
		I = LMi[i];

		for (j=0; j<M; ++j)
		{
			J = LMj[j];

			// only add values to lower-diagonal part of stiffness matrix
			if ((I>=J) && (J>=0)) 
			{
				pi = indices + pointers[J];
				l = pointers[J+1] - pointers[J];
				for (n=0; n<l; ++n) if (pi[n] == I) 
				{
					pd[pointers[J] + n] += ke[i][j];
					break;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
// FUNCTION : MatrixFactory::AssembleCompactUnSymm
// This function assembles the local stiffness matrix
// into the global stiffness matrix which is in compact column storage and
// the matrix is unsymmetric
//
void MatrixFactory::AssembleCompactUnSymm(CompactUnSymmMatrix& K, matrix& ke, vector<int>& LM)
{
	int i, j, I, J;

	const int N = ke.rows();

	for (i=0; i<N; ++i)
	{
		if ((I = LM[i])>=0)
		{
			for (j=0; j<N; ++j)
			{
				if ((J = LM[j]) >= 0) K.add(I,J, ke[i][j]);
			}
		}
	}
}

void MatrixFactory::AssembleCompactUnSymm(CompactUnSymmMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj)
{
	int i, j, I, J;

	const int N = ke.rows();
	const int M = ke.columns();

	for (i=0; i<N; ++i)
	{
		if ((I = LMi[i])>=0)
		{
			for (j=0; j<M; ++j)
			{
				if ((J = LMj[j]) >= 0) K.add(I,J, ke[i][j]);
			}
		}
	}
}
