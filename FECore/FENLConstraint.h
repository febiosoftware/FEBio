/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FESolver.h"
#include "FEStepComponent.h"
#include "FEGlobalVector.h"
#include "FEGlobalMatrix.h"
#include "FETimeInfo.h"
#include <vector>

//-----------------------------------------------------------------------------
// forward declaration of the model class
class FEModel;
class FELinearSystem;

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints enforced using an augmented Lagrangian method.

//! The constraint must provide a residual (force) contribution, its stiffness matrix,
//! and an augmentation function.
//!
class FECORE_API FENLConstraint : public FEStepComponent
{
	FECORE_SUPER_CLASS(FENLCONSTRAINT_ID)
	FECORE_BASE_CLASS(FENLConstraint);

public:
	FENLConstraint(FEModel* pfem);
	virtual ~FENLConstraint();

	// clone the constraint
	virtual void CopyFrom(FENLConstraint* plc) {}

	// initialize equations
	// Overridden by constraint classes that need to allocate more equations,
	// e.g. for Lagrange Multipliers.
	virtual int InitEquations(int neq) { return 0; }

public:
	// The LoadVector function evaluates the "forces" that contribute to the residual of the system
	virtual void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) = 0;

	// Evaluates the contriubtion to the stiffness matrix
	virtual void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) = 0;

	// Performs an augmentation step
	virtual bool Augment(int naug, const FETimeInfo& tp) { return true; }

	// Build the matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& M) = 0;

	// reset the state data
	virtual void Reset() {}

	// called at start of time step
	virtual void PrepStep() {}

	// update 
	using FEModelComponent::Update;
	virtual void Update(const std::vector<double>& ui) {}
	virtual void Update(const std::vector<double>& Ui, const std::vector<double>& ui) { Update(ui); }
	virtual void UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) {}
};
