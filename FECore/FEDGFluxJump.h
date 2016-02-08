#pragma once
#include "FESurfaceLoad.h"

//-----------------------------------------------------------------------------
// This class is a base class for evaluating jumps across inter-element boundaries.
class FEDGFluxJump : public FESurfaceLoad
{
public:
	//! constructor
	FEDGFluxJump(FEModel* pfem);

	//! one-time initialization
	bool Init();

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! calculate stiffness matrix
	void StiffnessMatrix(FESolver* psolver);
};
