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
#include "FEBiphasicSolver.h"
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FEBiphasicSolver to solve
// solute problems. 
class FEBIOMIX_API FEBiphasicSoluteSolver : public FEBiphasicSolver
{
public:
	//! con/descructor
	FEBiphasicSoluteSolver(FEModel* pfem);
	virtual ~FEBiphasicSoluteSolver(){}

	//! Initialize data structures
	bool Init() override;

	//! Initialize equations
	bool InitEquations() override;

	//! prepares the data for the first QN iteration
	void PrepStep() override;

	//! Performs a Newton-Raphson iteration
	bool Quasin() override;

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar) override;

public:
	//! Calculates residual (overridden from FEBiphasicSolver)
	bool Residual(vector<double>& R) override;

	//! calculates the global stiffness matrix (overridden from FESolidSolver2)
	bool StiffnessMatrix() override;

	//! update kinematics
	void UpdateKinematics(vector<double>& ui) override;

	//! Update solute data
	void UpdateSolute(vector<double>& ui);

protected:
	void GetConcentrationData(vector<double>& ci, vector<double>& ui, const int sol);

public:	// Parameters
	double	m_Ctol;			//!< concentration tolerance

public:
	// solute data
	vector< vector<double> >	m_ci;	//!< concentration increment vector
	vector< vector<double> >	m_Ci;	//!< Total concentration vector for iteration

	FEDofList	m_dofC;	//!< concentration dof
    FEDofList	m_dofD;	//!< shell concentration dof

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
