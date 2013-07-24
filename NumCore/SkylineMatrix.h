#pragma once
#include <FECore\SparseMatrix.h>

//=============================================================================
//! Implements a sparse matrix using the skyline storage

//! This class implements a symmetric sparse matrix where only the values
//! below the skyline are stored.

class SkylineMatrix : public SparseMatrix
{
public:
	SkylineMatrix();
	virtual ~SkylineMatrix();

	void Clear()
	{
		if (m_pd) delete [] m_pd; m_pd = 0;
		if (m_ppointers) delete [] m_ppointers; m_ppointers = 0;
	}

	void Create(SparseMatrixProfile& mp);

	void Assemble(matrix& ke, vector<int>& lm);

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj);

	void Create(double* pv, int* pp, int N);

	void add(int i, int j, double v)
	{
#ifdef _DEBUG
		if (j<i) { j ^= i; i ^= j; j ^= i; }
		int l = m_ppointers[j+1] - m_ppointers[j];
		assert(j-i<l);
#endif // _DEBUG

		// only add to the upper triangular part
		if (j>=i) m_pd[ m_ppointers[j] + j-i] += v;
	}

	void set(int i, int j, double v)
	{
		// only add to the upper triangular part
		if (j>=i) m_pd[ m_ppointers[j] + j-i] = v;
	}

	double get(int i, int j)
	{
		if (i>j) { i^= j; j ^= i; i ^= j; }

		int l = m_ppointers[j+1] - m_ppointers[j];
		if (j-i < l) return m_pd[ m_ppointers[j] + j-i];
		return 0;
	}

	double diag(int i)
	{
		return m_pd[ m_ppointers[i] ];
	}


	double* values() { return m_pd; }
	int* pointers() { return m_ppointers; }

protected:
	int*	m_ppointers;	// arrays of indices to diagonal elements
};
