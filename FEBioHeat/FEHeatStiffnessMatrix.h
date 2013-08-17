#pragma once
#include "FECore/FEGlobalMatrix.h"

//-----------------------------------------------------------------------------
//! The "stiffness" matrix for a heat-transfer problem
class FEHeatStiffnessMatrix : public FEGlobalMatrix
{
public:
	//! constructor
	FEHeatStiffnessMatrix(SparseMatrix* pK);

	//! destructor
	~FEHeatStiffnessMatrix();

	//! Construct the stiffness matrix from an FEModel
	bool Create(FEModel* pfem, int neq, bool breset);

protected:
	FEModel*	m_pfem;	//!< pointer to model
};
