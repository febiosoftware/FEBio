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
#include <FECore/FEGlobalVector.h>
#include <FECore/LinearSolver.h>
#include "FERigidSolver.h"
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEJacobiSolidSolver : public FESolver
{
public:
	//! constructor
	FEJacobiSolidSolver(FEModel* pfem);

	//! destructor
	virtual ~FEJacobiSolidSolver();

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
		bool PrepStep();

		//! Solve an analysis step
		bool SolveStep() override;

		//! solve equations
		void SolveEquations(std::vector<double>& u, std::vector<double>& R);
	//}

	//{ --- Stiffness matrix routines ---

		//! calculates the global stiffness matrix
		virtual bool StiffnessMatrix();

		//! contact stiffness
		void ContactStiffness(FELinearSystem& LS);

		//! calculates stiffness contributon of nonlinear constraints
		void NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp);
	//}

	//{ --- Residual routines ---

		//! Calculate the contact forces
		void ContactForces(FEGlobalVector& R);

		//! Calculates residual
		virtual bool Residual(vector<double>& R);

		//! Calculate nonlinear constraint forces
		void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

		//! Internal forces
		void InternalForces(FEGlobalVector& R);

		//! external forces
		void ExternalForces(FEGlobalVector& R);
	//}

private:
	bool CalculateMassMatrix();

	bool AllocateLinearSystem();

	bool DoAugmentations();

	//! reform the stiffness matrix
	bool ReformStiffness();

	//! recalculates the shape of the stiffness matrix
	bool CreateStiffness(bool breset);

public:
	// solver parameters
	double	m_Etol;			//!< energy tolerance
	double	m_Rtol;			//!< residual tolerance
	double	m_Dtol;			//!< displacement tolerance
	double	m_Rmin;			//!< min residual value
	double	m_Rmax;			//!< max residual value
	int		m_maxref;		//!< max nr of reformations per time step
	bool	m_bdivreform;		//!< reform when diverging
	bool	m_breformAugment;	//!< reform after each (failed) augmentations
	bool	m_breformtimestep;	//!< reform at start of time step
	bool	m_bforceReform;		//!< forces a reform in QNInit
	bool	m_bdoreforms;		//!< do reformations

	// equation numbers
	int		m_nreq;			//!< start of rigid body equations

public:
	vector<double> m_Fr;	//!< nodal reaction forces
	vector<double> m_Fint;	//!< internal load vector
	vector<double> m_Fext;	//!< external load vector
	vector<double> m_R0;	//!< residual at iteration i-1
	vector<double> m_R1;	//!< residual at iteration i
	vector<double> m_ui;	//!< solution increment vector
	vector<double> m_Ut;	//!< total solution vector
	vector<double> m_Ui;	//!< total solution increments of current time step
	vector<double> m_up;	//!< solution increment of previous iteration
	vector<double> m_Fd;	//!< residual correction due to prescribed degrees of freedom
	vector<double> m_Mi;	//!< inverse of lumped mass matrix

    // generalized alpha method (for dynamic analyses)
    double  m_rhoi;         //!< spectral radius
    double  m_alphaf;       //!< alpha step for Y={v,e}
    double  m_alpham;       //!< alpha step for Ydot={dv/dt,de/dt}
	double	m_alpha;		//!< Newmark parameter alpha (force integration)
	double	m_beta;			//!< Newmark parameter beta (displacement integration)
	double	m_gamma;		//!< Newmark parameter gamme (velocity integration)

private:
	LinearSolver*	m_plinsolve;	//!< the linear solver (NOTE: WE don't really need it, but we use it to allocate a sparse matrix)
	FEGlobalMatrix* m_pK;			//!< global stiffness matrix
	bool			m_breshape;		//!< Matrix reshape flag

	// Error handling
	bool	m_bzero_diagonal;	//!< check for zero diagonals
	double	m_zero_tol;			//!< tolerance for zero diagonal

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
