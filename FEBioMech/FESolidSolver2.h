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
//! The FESolidSolver2 class solves large deformation solid mechanics problems
//! It can deal with quasi-static and dynamic problems
//! 
class FEBIOMECH_API FESolidSolver2 : public FENewtonSolver
{
	enum ARC_LENGTH_METHOD {
		NONE,
		CRISFIELD,
	};

public:
	//! constructor
	FESolidSolver2(FEModel* pfem);

	//! destructor
	virtual ~FESolidSolver2();

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar) override;

	//! Initializes data structures
	bool Init() override;

	//! initialize the step
	bool InitStep(double time) override;

	//! Initialize linear equation system
	bool InitEquations() override;
	bool InitEquations2() override;

    //! Generate warnings if needed
    void SolverWarnings();

	//! Return the rigid solver
	FERigidSolver* GetRigidSolver();

public:
	//{ --- evaluation and update ---
		//! Perform an update
		void Update(vector<double>& ui) override;

		//! perform an updated where ui also contains displacement increments of prescribed displacements
		//! NOTE: This is a temporary hack that is only by the JFNKMatrix
		void Update2(const vector<double>& ui) override;

		//! update nodal positions, velocities, accelerations, etc.
		virtual void UpdateKinematics(vector<double>& ui);

		//! Update EAS
		void UpdateEAS(vector<double>& ui);
		void UpdateIncrementsEAS(vector<double>& ui, const bool binc);

		//! update DOF increments
		virtual void UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap);
	//}

	//{ --- Solution functions ---

		//! prepares the data for the first QN iteration
		void PrepStep() override;

		//! Performs a Newton-Raphson iteration
		bool Quasin() override;

		//! Apply arc-length
		void DoArcLength();
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

		//! Calculate the contact forces
		void ContactForces(FEGlobalVector& R);

		//! Calculates residual
		virtual bool Residual(vector<double>& R) override;

		//! Calculate nonlinear constraint forces
		void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

		//! Internal forces
		void InternalForces(FEGlobalVector& R);

		//! external forces
		void ExternalForces(FEGlobalVector& R);
	//}

public:
	// convergence tolerances
	double	m_Dtol;			//!< displacement tolerance

	bool	m_logSolve;		//!< flag to use Aggarwal's log method

	// equation numbers
	int		m_nreq;			//!< start of rigid body equations

public:
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_Fint;	//!< internal load vector
	vector<double> m_Fext;	//!< external load vector
	vector<double> m_Uip;	//!< previous converged displacement increment

    // generalized alpha method (for dynamic analyses)
    double  m_rhoi;         //!< spectral radius
    double  m_alphaf;       //!< alpha step for Y={v,e}
    double  m_alpham;       //!< alpha step for Ydot={∂v/∂t,∂e/∂t}
	double	m_alpha;		//!< Newmark parameter alpha (force integration)
	double	m_beta;			//!< Newmark parameter beta (displacement integration)
	double	m_gamma;		//!< Newmark parameter gamme (velocity integration)

	// arc-length parameters
	int		m_arcLength;	//!< arc-length method flag (0 = off, 1 = Crisfield)
	double	m_al_scale;		//!< arc-length scaling parameter (i.e. psi).
	double	m_al_lam;		//!< current arc-length lambda
	double	m_al_inc;		//!< arc-length lambda increment at current timestep
	double	m_al_ds;		//!< arc-length constraint
	double	m_al_gamma;		//!< acr-length increment at current iteration

protected:
	FEDofList	m_dofU, m_dofV;
	FEDofList	m_dofSQ;
	FEDofList	m_dofRQ;
	FEDofList	m_dofSU, m_dofSV, m_dofSA;
    
protected:
    FERigidSolverNew	m_rigidSolver;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
