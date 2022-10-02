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
#include "matrix.h"
#include "FENewtonStrategy.h"

//-----------------------------------------------------------------------------
//! This class implements the Broyden quasi-newton strategy. 
class FECORE_API FEBroydenStrategy : public FENewtonStrategy
{
public:
	//! constructor
	FEBroydenStrategy(FEModel* fem);

	//! Initialization
	bool Init() override;

	//! perform a quasi-Newton udpate
	bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1) override;

	//! solve the equations
	void SolveEquations(vector<double>& x, vector<double>& b) override;

	//! Presolve update
	virtual void PreSolveUpdate() override;

private:
	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;	//!< pointer to linear solver
	int				m_neq;			//!< number of equations

	bool		m_bnewStep;

	// Broyden update vectors
	matrix			m_R;		//!< Broyden update vector "r"
	matrix			m_D;		//!< Broydeb update vector "delta"
	vector<double>	m_rho;		//!< temp vectors for calculating Broyden update vectors
	vector<double>	m_q;		//!< temp storage for q

	DECLARE_FECORE_CLASS();
};
