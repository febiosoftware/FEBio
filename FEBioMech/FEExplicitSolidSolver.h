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
#include "FECore/FESolver.h"
#include "FECore/FEGlobalVector.h"
#include <FECore/FETimeInfo.h>
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
//! This class implements a nonlinear explicit solver for solid mechanics
//! problems.
class FEExplicitSolidSolver : public FESolver
{
	enum MassLumpingMethod
	{
		NO_MASS_LUMPING,	// use consistent mass matrix
		ROW_SUM_LUMPING,	// use simple row-sum lumping
		HRZ_LUMPING			// use Hinton-Rock-Zienkiewicz lumping
	};

public:
	//! constructor
	FEExplicitSolidSolver(FEModel* pfem);

	//! destructor
	virtual ~FEExplicitSolidSolver() {}

public:
	//! Data initialization
	bool Init() override;

	//! clean up
	void Clean() override;

	//! Solve an analysis step
	bool SolveStep() override;

	//! Update data
	void Update(vector<double>& ui) override;

	//! Serialize data
	void Serialize(DumpStream& ar) override;

public:
	//! update kinematics
	void UpdateKinematics(vector<double>& ui);

	//! Update rigid bodies 
	void UpdateRigidBodies(vector<double>& ui);

	//! solve the step
	bool DoSolve();

	void PrepStep();

	bool Residual(vector<double>& R);

	void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

	void ContactForces(FEGlobalVector& R);

private:
	bool CalculateMassMatrix();

public:
	int			m_mass_lumping;	//!< specify mass lumping method
	double		m_dyn_damping;	//!< velocity damping for the explicit solver

public:
	// equation numbers
	int		m_nreq;			//!< start of rigid body equations

	vector<double> m_Mi;	//!< inverse mass vector for explicit analysis
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_Ut;	//!< Total dispalcement vector at time t (incl all previous timesteps)

	vector<double> m_ui;	//!< displacement increment vector

	vector<double> m_R0;	//!< residual at iteration i-1
	vector<double> m_R1;	//!< residual at iteration i

protected:
	FEDofList	m_dofU, m_dofV, m_dofSQ, m_dofRQ;
	FEDofList	m_dofSU, m_dofSV, m_dofSA;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
