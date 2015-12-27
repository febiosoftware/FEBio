#pragma once
#include "FESolver.h"
#include "DumpFile.h"
#include "FEModelComponent.h"
#include "FEGlobalVector.h"
#include "FEGlobalMatrix.h"
#include "FESurface.h"
#include "FETypes.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// forward declaration of the model class
class FEModel;

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints enforced using an augmented Lagrangian method.

//! The constraint must provide a residual (force) contribution, its stiffness matrix
//! and an augmentation function.
//!
class FENLConstraint : public FEModelComponent
{
public:
	FENLConstraint(FEModel* pfem);
	virtual ~FENLConstraint();

public:
	virtual void Residual(FEGlobalVector& R, const FETimePoint& tp) = 0;
	virtual void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp) = 0;
	virtual bool Augment(int naug, const FETimePoint& tp) = 0;
	virtual void Serialize(DumpFile& ar) = 0;
	virtual void ShallowCopy(DumpStream& dmp, bool bsave) = 0;
	virtual void CopyFrom(FENLConstraint* plc) {}
	virtual void BuildMatrixProfile(FEGlobalMatrix& M) = 0;

	// update state
	virtual void Reset() {}
	virtual void Update(const FETimePoint& tp) {}

	virtual FESurface* GetSurface(const char* sz) { return 0; }
};
