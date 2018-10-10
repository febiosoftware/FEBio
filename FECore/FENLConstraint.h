#pragma once
#include "FESolver.h"
#include "FEModelComponent.h"
#include "FEGlobalVector.h"
#include "FEGlobalMatrix.h"
#include "FETimeInfo.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// forward declaration of the model class
class FEModel;

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints enforced using an augmented Lagrangian method.

//! The constraint must provide a residual (force) contribution, its stiffness matrix,
//! and an augmentation function.
//!
class FECORE_API FENLConstraint : public FEModelComponent
{
	DECLARE_SUPER_CLASS(FENLCONSTRAINT_ID);

public:
	FENLConstraint(FEModel* pfem);
	virtual ~FENLConstraint();

	// clone the constraint
	virtual void CopyFrom(FENLConstraint* plc) {}

public:
	// The Residual function evaluates the "forces" that contribute to the residual of the system
	virtual void Residual(FEGlobalVector& R, const FETimeInfo& tp) = 0;

	// Evaluates the contriubtion to the stiffness matrix
	virtual void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) = 0;

	// Performs an augmentation step
	virtual bool Augment(int naug, const FETimeInfo& tp) = 0;

	// Build the matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& M) = 0;

	// Update state
	virtual void Update(int niter, const FETimeInfo& tp) {}

	// reset the state data
	virtual void Reset() {}
};
