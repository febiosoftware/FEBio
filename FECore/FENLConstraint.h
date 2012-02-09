#pragma once
#include "FENLSolver.h"
#include "DumpFile.h"
#include "FEParameterList.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// types of nonlinear constraints
enum FENLConstraint_Type {
	FE_POINT_CONSTRAINT,
	FE_LINEAR_CONSTRAINT,
	FE_RIGID_JOINT
};

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
	FENLConstraint(FEModel* pfem, int ntype) : m_pfem(pfem), m_ntype(ntype) {}
	virtual ~FENLConstraint(){}

	int Type() const { return m_ntype; }

public:
	virtual void Init() = 0;
	virtual void Residual(FENLSolver* psolver, vector<double>& R) = 0;
	virtual void StiffnessMatrix(FENLSolver* psolver) = 0;
	virtual bool Augment(int naug) = 0;
	virtual void Serialize(DumpFile& ar) = 0;

	// update state
	virtual void Update() {}

protected:
	FEModel*	m_pfem;
	int			m_ntype;
};
