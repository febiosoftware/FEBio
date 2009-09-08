// MatrixFactory.h: interface for the MatrixFactory class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIXFACTORY_H__9E155491_5839_40CE_B59D_1BC0C9383962__INCLUDED_)
#define AFX_MATRIXFACTORY_H__9E155491_5839_40CE_B59D_1BC0C9383962__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "SparseMatrix.h"
#include "vector.h"
#include "matrix.h"
#include "MatrixProfile.h"

////////////////////////////////////////////////////////////////////////////////////
// CLASS : MatrixFactory
// The matrix factory is the link between the FE code and the different matrix
// classes. Since the FE code has no idea what type of matrix it's using we need
// to provide some mechanism to make sure data is placed correctly in the 
// stiffness matrix. The MatrixFactory provides the necessary facilities to make
// this happen.
// In the first place the matrix factory is used when initial allocation for the
// stiffness matrix is requested. The CreateMatrix function will take care of this.
//

class MatrixFactory  
{
public:
	// constructor/destructor
	MatrixFactory();
	virtual ~MatrixFactory();

	// Allocates storage for sparse matrix A
	bool CreateMatrix(SparseMatrix* pA, vector< vector<int> >& LM, int neq);

	// creates a matrix from a MatrixProfile
	bool CreateMatrix(SparseMatrix* pA, MatrixProfile& mp);

	// Assembles element stiffness matrix into global stiffness matrix
	void Assemble(SparseMatrix& K, matrix& ke, vector<int>& LM);

	// More general version of assembly routine
	void Assemble(SparseMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj);

private:
	// creates a full matrix
	bool CreateFull(FullMatrix* pA, vector< vector<int> >& LM, int neq);
	bool CreateFull(FullMatrix* pA, MatrixProfile& mp);

	// creates a symmetric skyline matrix
	bool CreateSkyline(SkylineMatrix* pA, vector< vector<int> >& LM, int neq);
	bool CreateSkyline(SkylineMatrix* pA, MatrixProfile& mp);

	// creates a symmetric compact (Harwell-Boeing) matrix
	bool CreateCompact(CompactSymmMatrix* pA, vector< vector<int> >& LM, int neq);
	bool CreateCompact(CompactSymmMatrix* pA, MatrixProfile& mp);

	// creates a compact unsymmetric (Harwell-Boeing) matrix
	bool CreateCompactUnSymm(CompactUnSymmMatrix* pA, vector< vector<int> >& LM, int neq);
	bool CreateCompactUnSymm(CompactUnSymmMatrix* pA, MatrixProfile& mp);

	// Assembler routine for full matrices
	void AssembleFull(FullMatrix& K, matrix& ke, vector<int>& LM);
	void AssembleFull(FullMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj);

	// Assembler routine for skyline matrices
	void AssembleSkyline(SkylineMatrix& K, matrix& ke, vector<int>& LM);
	void AssembleSkyline(SkylineMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj);

	// Assembler routine for compact matrices
	void AssembleCompact(CompactSymmMatrix& K, matrix& ke, vector<int>& LM);
	void AssembleCompact(CompactSymmMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj);

	// Assembler routine for unsymmetric compact matrices
	void AssembleCompactUnSymm(CompactUnSymmMatrix& K, matrix& ke, vector<int>& LM);
	void AssembleCompactUnSymm(CompactUnSymmMatrix& K, matrix& ke, vector<int>& LMi, vector<int>& LMj);
};

#endif // !defined(AFX_MATRIXFACTORY_H__9E155491_5839_40CE_B59D_1BC0C9383962__INCLUDED_)
