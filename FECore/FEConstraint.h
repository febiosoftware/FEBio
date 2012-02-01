#pragma once
#include "FENLSolver.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// forward declaration of the model class
class FEModel;

//-----------------------------------------------------------------------------
// Base class for nonlinear constraints enforced using an augmented Lagrangian method.
// TODO: change the name to FENLConstraint
class FEConstraint
{
public:
	FEConstraint(FEModel* pfem) : m_pfem(pfem) {}
	virtual ~FEConstraint(){}

public:
	virtual void Init() = 0;
	virtual void Residual(FENLSolver* psolver, vector<double>& R) = 0;
	virtual void StiffnessMatrix(FENLSolver* psolver) = 0;
	virtual bool Augment(int naug) = 0;

protected:
	FEModel*	m_pfem;
};
