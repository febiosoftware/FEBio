#pragma once

#include "SparseMatrix.h"
#include "FENLSolver.h"
#include <vector>

//-----------------------------------------------------------------------------
//! This class implements a global system matrix.

//! The global system matrix is usually created by the discretization of the FE 
//! equations into a linear system of equations. The structure of it depends greatly
//! on the element connectivity and usually results in a sparse matrix structure. 
//! FEBio supports several sparse matrix structure (Compact, Skyline, etc.) and to 
//! simplify the creation of the specific matrix structure, the FEGlobalMatrix offers
//! functionality to create the global matrix structure without the need to know 
//! what particular sparse matrix format is used by the linear solver.

//! \todo I think the SparseMatrixProfile can handle all of the build functions.

class FEGlobalMatrix
{
protected:
	enum { MAX_LM_SIZE = 4096 };

public:
	//! constructor
	FEGlobalMatrix(SparseMatrix* pK);

	//! destructor
	virtual ~FEGlobalMatrix();

	//! construct the stiffness matrix from a FEM object
	virtual bool Create(FENLSolver* psolver, int neq, bool breset) = 0;

	//! clears the sparse matrix that stores the stiffness matrix
	void Clear();

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

	//! return a pointer to the sparse matrix
	SparseMatrix* GetSparseMatrixPtr() { return m_pA; }

	//! zero the sparse matrix
	void Zero() { m_pA->zero(); }

protected:
	void build_begin(int neq);
	void build_add(std::vector<int>& lm);
	void build_end();
	void build_flush();

protected:
	SparseMatrix*	m_pA;	//!< the actual global stiffness matrix

	// The following data structures are used to incrementally
	// build the profile of the sparse matrix
	SparseMatrixProfile*	m_pMP;		//!< profile of sparse matrix
	vector< vector<int> >			m_LM;		//!< used for building the stiffness matrix
	int	m_nlm;									//!< nr of elements in m_LM array
};
