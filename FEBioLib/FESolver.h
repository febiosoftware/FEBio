// FESolver.h: interface for the FESolidSolver class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_)
#define AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEBioLib/FEStiffnessMatrix.h"
#include "FEBioLib/Timer.h"
#include "FECore/FEException.h"
#include "FECore/FENLSolver.h"
using namespace NumCore;

//-----------------------------------------------------------------------------
//! The solver class.

//! This class is responsible for solving a single time step in the FE
//! analysis, possibly using results from previous solutions.
//! In the future this class might become a base-class for different FE solvers.
//! For instances there could be a different solver for quasi-static, dynamic and
//! eigenvalue problems.
class FESolver : public FENLSolver
{
public:
	FESolver(FEModel& fem);
	virtual ~FESolver();

	virtual bool Init();
	virtual void Clean();

public:
	//! recalculates the shape of the stiffness matrix
	bool CreateStiffness(bool breset);

	//! return pointer to stiffness matrix
	FEStiffnessMatrix* GetStiffnessMatrix() { return m_pK; }

public:
	//! assemble into global residual (TODO: this is only used by rigid joints)
//	void AssembleResidual(vector<int>& lm, vector<double>& fe, vector<double>& R);

	//! assemble global stiffness matrix (TODO: this is only used by rigid joints)
	void AssembleStiffness(vector<int>& elm, matrix& ke);

public:
	// timers
	Timer	m_SolverTime;	//!< time spent in solver

	// linear solver
	LinearSolver*		m_plinsolve;	//!< the linear solver

	// global stiffness matrix
	FEStiffnessMatrix*	m_pK;		//!< global stiffness matrix
	int					m_neq;		//!< number of equations
};

#endif // !defined(AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_)
