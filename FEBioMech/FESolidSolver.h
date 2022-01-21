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

#include "FECore/FENewtonSolver.h"
#include <FECore/FETimeInfo.h>
#include "FECore/FEGlobalVector.h"
#include "FERigidSolver.h"
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
//! The FESolidSolver class solves large deformation solid mechanics problems
//! It can deal with quasi-static and dynamic problems
//! 
class FESolidSolver : public FENewtonSolver
{
public:
	//! constructor
	FESolidSolver(FEModel* pfem);

	//! destructor
	virtual ~FESolidSolver();

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar) override;

	//! Initializes data structures
	bool Init() override;

	//! Initialize linear equation system
	bool InitEquations() override;

	// get the rigid solver
	FERigidSolver* GetRigidSolver() { return &m_rigidSolver; }

public:
	//{ --- evaluation and update ---
		//! Perform an update
		void Update(vector<double>& ui) override;
	//}

	//{ --- Solution functions ---

		//! prepares the data for the first QN iteration
		void PrepStep() override;

		//! Performs a Newton-Raphson iteration
		bool Quasin() override;

		//! update nodal positions, velocities, accelerations, etc.
		virtual void UpdateKinematics(vector<double>& ui);
	//}

	//{ --- Stiffness matrix routines ---

		//! calculates the global stiffness matrix
		virtual bool StiffnessMatrix() override;

		//! contact stiffness
		void ContactStiffness(FELinearSystem& LS);

		//! calculates stiffness contributon of nonlinear constraints
		void NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp);
	//}

	//{ --- Residual routines ---

		//! Calculate inertial forces for dynamic problems
		void InertialForces(FEGlobalVector& R);

		//! Calculate the contact forces
		void ContactForces(FEGlobalVector& R);

		//! Calculates residual
		virtual bool Residual(vector<double>& R) override;

		//! Calculate nonlinear constraint forces
		void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);
	//}

public:
	// convergence tolerances
	double	m_Dtol;			//!< displacement tolerance

	// equation numbers
	int		m_nreq;			//!< start of rigid body equations

	// Newmark parameters (for dynamic analyses)
	double	m_beta;			//!< Newmark parameter beta (displacement integration)
	double	m_gamma;		//!< Newmark parameter gamme (velocity integration)

public:
	vector<double> m_Fr;	//!< nodal reaction forces

public:
	bool	m_bnew_update;	//!< use new rigid body update algorithm

protected:
	FEDofList	m_dofU, m_dofV, m_dofSQ, m_dofRQ;

protected:
	FERigidSolverOld	m_rigidSolver;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
