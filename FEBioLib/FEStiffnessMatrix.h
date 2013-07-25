// FEStiffnessMatrix.h: interface for the FEStiffnessMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESTIFFNESSMATRIX_H__8E7BEF6B_A12D_4C74_9C88_3ADE0141B981__INCLUDED_)
#define AFX_FESTIFFNESSMATRIX_H__8E7BEF6B_A12D_4C74_9C88_3ADE0141B981__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <FECore\FEGlobalMatrix.h>
#include <FECore\FEModel.h>

//-----------------------------------------------------------------------------
//! The FEStiffnessmatrix class stores the global stiffness matrix of the FE
//! problem. It also contains the algorithm to construct the stiffness matrix
//! from a FEM object. The actual matrix data is stored in a SparseMatrix class.

class FEStiffnessMatrix : public FEGlobalMatrix
{
public:
	//! constructor
	FEStiffnessMatrix(SparseMatrix* pK);

	//! destructor
	virtual ~FEStiffnessMatrix();

	//! construct the stiffness matrix from a FEM object
	bool Create(FEModel* pfem, int neq, bool breset);

protected:
	void AddContactInterface(FEContactInterface* pci);

protected:
	FEModel*		m_pfem;	//!< pointer to model
};

#endif // !defined(AFX_FESTIFFNESSMATRIX_H__8E7BEF6B_A12D_4C74_9C88_3ADE0141B981__INCLUDED_)
