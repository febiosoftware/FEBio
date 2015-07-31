#pragma once
#include "FEModelComponent.h"
#include "FEGlobalVector.h"
#include "FETypes.h"

//-----------------------------------------------------------------------------
class FESolver;

//-----------------------------------------------------------------------------
//! This class is the base class for all classes that affect the state of the model
//! and contribute directly to the residual and the global stiffness matrix. This
//! includes most boundary loads, body loads, contact, etc.
class FEModelLoad : public FEModelComponent
{
public:
	//! constructor
	FEModelLoad(SUPER_CLASS_ID sid, FEModel* pfem);

public:
	// all classes derived from this base class must implement
	// the following functions.

	//! evaluate the contribution to the residual
	virtual void Residual(FEGlobalVector& R, const FETimePoint& tp) = 0;

	//! evaluate the contribution to the global stiffness matrix
	virtual void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp) = 0;
};
