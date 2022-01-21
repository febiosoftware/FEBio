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
#include <FEBioMech/FESolidSolver2.h>
#include <FECore/FEElementTraits.h>
#include <FECore/FEDofList.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver to solve
// biphasic problems. 
class FEBIOMIX_API FEBiphasicSolver : public FESolidSolver2
{
public:
	//! constructor
	FEBiphasicSolver(FEModel* pfem);
	virtual ~FEBiphasicSolver();

	//! Initialize data structures
	bool Init() override;

	//! Initialize linear equation system
	bool InitEquations() override;

	//! prepares the data for the first QN iteration
	void PrepStep() override;

	//! Performs a Newton-Raphson iteration
	bool Quasin() override;

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar) override;

public:
	//! update contact
	void UpdateModel() override;

	//! update kinematics
	void UpdateKinematics(vector<double>& ui) override;

	//! Update poroelastic data
	void UpdatePoro(vector<double>& ui);

public:

	//! Calculates residual (overridden from FESolidSolver2)
	bool Residual(vector<double>& R) override;

	//! calculates the global stiffness matrix (overridden from FESolidSolver2)
	bool StiffnessMatrix() override;

protected:
	void GetDisplacementData(vector<double>& di, vector<double>& ui);
	void GetPressureData(vector<double>& pi, vector<double>& ui);

public:
	// additional convergence norms
	double	m_Ptol;			//!< pressure tolerance

	// biphasic formulation
	int		m_biphasicFormulation;	// = 0: standard, =1: mixed (linear pressure)

	// equation numbers
	int		m_ndeq;				//!< number of equations related to displacement dofs
	int		m_npeq;				//!< number of equations related to pressure dofs
	vector<int>		m_nceq;	//!< number of equations related to concentration dofs

	// poro data
	vector<double>	m_di;	//!< displacement increment vector
	vector<double>	m_Di;	//!< total displacement vector for iteration
	vector<double>	m_pi;	//!< pressure increment vector
	vector<double>	m_Pi;	//!< Total pressure vector for iteration

protected:
	FEDofList	m_dofP;	//!< pressure dof index
	FEDofList	m_dofQ;	//!< shell pressure dof index
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
