#pragma once
#include <vector>
#include "SparseMatrix.h"

namespace NumCore {

//----------------------------------------------------------------------------
//! The NonLinearSystem describes the interface for all non-linear systems.

//! Each non-linear system has four responsibilities:
//! 1. Evaluate its current state
//! 2. Calculate its current jacobian
//! 3. Update its current state given a solution increment
//! 4. Decide whether the current state is considered converged
class NonLinearSystem
{
public:
	//! Class destructor
	virtual ~NonLinearSystem(void);

	//! overide function to evaluate current state
	virtual void Evaluate(std::vector<double>& F) = 0;

	//! override function to calculate jacobian matrix
	virtual void Jacobian(SparseMatrix& K) = 0;

	//! override function to update state of system
	virtual void Update(std::vector<double>& u) = 0;

	//! override to decide whether solution is converged
	virtual bool Converged() = 0;
};

} // namespace NumCore
