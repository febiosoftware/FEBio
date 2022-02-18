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
#include "FELinearConstraint.h"
#include "table.h"

class FEGlobalMatrix;
class matrix;

//-----------------------------------------------------------------------------
// This class helps manage all the linear constraints
class FECORE_API FELinearConstraintManager
{
public:
	FELinearConstraintManager(FEModel* fem);
	~FELinearConstraintManager();

	// Clear all constraints
	void Clear();

	// copy data
	void CopyFrom(const FELinearConstraintManager& lcm);

	// serialize linear constraints
	void Serialize(DumpStream& ar);

	// build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& G);

	// add the linear constraint
	void AddLinearConstraint(FELinearConstraint* lc);

	// return number of linear constraints
	int LinearConstraints() const;

	// return linear constraint
	const FELinearConstraint& LinearConstraint(int i) const;

	// return linear constraint
	FELinearConstraint& LinearConstraint(int i);

	//! remove a linear constraint
	void RemoveLinearConstraint(int i);

public:
	// one-time initialization
	bool Initialize();

	// activation
	bool Activate();

	// assemble element residual into global residual
	void AssembleResidual(vector<double>& R, vector<int>& en, vector<int>& elm, vector<double>& fe);

	// assemble element matrix into (reduced) global matrix
	void AssembleStiffness(FEGlobalMatrix& K, vector<double>& R, vector<double>& ui, const vector<int>& en, const vector<int>& lmi, const vector<int>& lmj, const matrix& ke);

	// called before the first reformation for each time step
	void PrepStep();

	// update nodal variables
	void Update();

protected:
	void InitTable();

private:
	FEModel* m_fem;
	vector<FELinearConstraint*>	m_LinC;		//!< linear constraints data
	table<int>					m_LCT;		//!< linear constraint table
	vector<double>				m_up;		//!< the inhomogenous component of the linear constraint
};
