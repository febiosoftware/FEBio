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
#include <FECore/FENewtonSolver.h>
#include <FECore/FEElementTraits.h>
#include <FECore/FEDofList.h>
#include <FEBioMech/FERigidSolver.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
// This class adds additional functionality to the FESolidSolver to solve
// biphasic problems. 
class FEBIOMIX_API FEBiphasicSolver : public FENewtonSolver
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

	//! Generate warnings if needed
	void SolverWarnings();

public:
	void Update(vector<double>& ui) override;

	//! update contact
	void UpdateModel() override;

	//! update kinematics
	virtual void UpdateKinematics(vector<double>& ui);

	//! Update EAS
	void UpdateEAS(vector<double>& ui);
	void UpdateIncrementsEAS(vector<double>& ui, const bool binc);

	//! Update poroelastic data
	void UpdatePoro(vector<double>& ui);

public:

	//! Calculates residual (overridden from FESolidSolver2)
	bool Residual(vector<double>& R) override;

	//! calculates the global stiffness matrix (overridden from FESolidSolver2)
	bool StiffnessMatrix() override;

    //! Internal forces
    void InternalForces(FEGlobalVector& R);
    
    //! external forces
    void ExternalForces(FEGlobalVector& R);

	void ContactForces(FEGlobalVector& R);

	void ContactStiffness(FELinearSystem& LS);

	void NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp);

	void NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp);

protected:
	void GetDisplacementData(vector<double>& di, vector<double>& ui);
	void GetPressureData(vector<double>& pi, vector<double>& ui);

public:
	// additional convergence norms
	double	m_Dtol;			//!< displacement tolerance
	double	m_Ptol;			//!< pressure tolerance
    double  m_Ctol;         //!< needed only for biphasic-solute analyses

	// biphasic formulation
	int		m_biphasicFormulation;	// = 0: standard, =1: mixed (linear pressure)

	// equation numbers
	int		m_ndeq;				//!< number of equations related to displacement dofs
	int		m_npeq;				//!< number of equations related to pressure dofs
	int		m_nreq;			//!< start of rigid body equations
	vector<int>		m_nceq;	//!< number of equations related to concentration dofs

	vector<double> m_Fr;	//!< nodal reaction forces

	// poro data
	vector<double>	m_di;	//!< displacement increment vector
	vector<double>	m_Di;	//!< total displacement vector for iteration
	vector<double>	m_pi;	//!< pressure increment vector
	vector<double>	m_Pi;	//!< Total pressure vector for iteration

protected:
	FEDofList	m_dofU, m_dofV;
	FEDofList	m_dofRQ;
	FEDofList	m_dofSU, m_dofSV, m_dofSA;
	FEDofList	m_dofP;	//!< pressure dof index
	FEDofList	m_dofSP;	//!< shell pressure dof index

	// obsolete parameters
	double  m_rhoi = -2;
	double	m_alpha = 1;
	double	m_beta = 0.25;
	double	m_gamma = 0.5;
	bool	m_logSolve = false;
	int		m_arcLength = 0;
	double	m_al_scale = 0;

protected:
	FERigidSolverNew	m_rigidSolver;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
