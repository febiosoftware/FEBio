#pragma once
#include "FENLSolver.h"
#include "DumpFile.h"
#include "FEParameterList.h"
#include "FEGlobalVector.h"
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
class FENLConstraint : public FEParamContainer
{
public:
	FENLConstraint(FEModel* pfem) : m_pfem(pfem) {}
	virtual ~FENLConstraint(){}

public:
	virtual void Init() = 0;
	virtual void Residual(FEGlobalVector& R) = 0;
	virtual void StiffnessMatrix(FENLSolver* psolver) = 0;
	virtual bool Augment(int naug) = 0;
	virtual void Serialize(DumpFile& ar) = 0;

	// update state
	virtual void Update() {}

protected:
	FEModel*	m_pfem;
};
