// FEStiffnessMatrix.h: interface for the FEStiffnessMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESTIFFNESSMATRIX_H__8E7BEF6B_A12D_4C74_9C88_3ADE0141B981__INCLUDED_)
#define AFX_FESTIFFNESSMATRIX_H__8E7BEF6B_A12D_4C74_9C88_3ADE0141B981__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/SparseMatrix.h"
using namespace FECore;

class FEM;
class FESolver;

//-----------------------------------------------------------------------------
//! The FEStiffnessmatrix class stores the global stiffness matrix of the FE
//! problem. It also contains the algorithm to construct the stiffness matrix
//! from a FEM object. The actual matrix data is stored in a SparseMatrix class.

class FEStiffnessMatrix  
{
public:
	//! constructor
	FEStiffnessMatrix(SparseMatrix* pK);

	//! destructor
	virtual ~FEStiffnessMatrix();

	//! clears the sparse matrix that stores the stiffness matrix
	void Clear() { if (m_pA) m_pA->Clear(); }

	//! construct the stiffness matrix from a FEM object
	bool Create(FESolver* psolver, int neq, bool breset);

	//! assemble an element stiffness matrix into the global stiffness matrix
	void Assemble(matrix& ke, vector<int>& lm) { m_pA->Assemble(ke, lm); }

	//! more general assembly routine
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj) { m_pA->Assemble(ke, lmi, lmj); }

	//! return the nonzeroes in the sparse matrix
	int NonZeroes() { return m_pA->NonZeroes(); }

	//! return the number of rows
	int Rows() { return m_pA->Size(); }

	//! converts a FEStiffnessMatrix to a SparseMatrix
	operator SparseMatrix* () { return m_pA; }

	//! converts a FEStiffnessMatrix to a SparseMatrix
	operator SparseMatrix& () { return *m_pA;}

	SparseMatrix* GetSparseMatrixPtr() { return m_pA; }

	void Zero() { m_pA->zero(); }

protected:
	void build_begin(int neq);
	void build_add(vector<int>& lm);
	void build_end();
	void build_flush();

protected:
	enum { MAX_LM_SIZE = 4096 };

protected:
	SparseMatrix*	m_pA;	//!< the actual global stiffness matrix

	// The following data structures are used to incrementally
	// build the profile of the sparse matrix
	FEM*	m_pfem;					//!< pointer to the FEM object
	SparseMatrixProfile*	m_pMP;			//!< profile of sparse matrix
	vector< vector<int> >	m_LM;	//!< used for building the stiffness matrix
	int	m_nlm;						//!< nr of elements in m_LM array
};

#endif // !defined(AFX_FESTIFFNESSMATRIX_H__8E7BEF6B_A12D_4C74_9C88_3ADE0141B981__INCLUDED_)
