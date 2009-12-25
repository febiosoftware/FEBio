// FESolver.h: interface for the FESolidSolver class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_)
#define AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/LinearSolver.h"
#include "FEStiffnessMatrix.h"
#include "fem.h"
#include "Timer.h"
#include "FEException.h"
#include "Interrupt.h"
#include "BFGSSolver.h"

//-----------------------------------------------------------------------------
//! The solver class.

//! This class is responsible for solving a single time step in the FE
//! analysis, possibly using results from previous solutions.
//! In the future this class might become a base-class for different FE solvers.
//! For instances there could be a different solver for quasi-static, dynamic and
//! eigenvalue problems.
//! This class is derived from Interruptable, because we want to be able to capture
//! the CTRL+C interruptions

class FESolver : public Interruptable
{
public:
	FESolver(FEM& fem);
	virtual ~FESolver();

	virtual bool Init();
	virtual bool SolveStep(double time) = 0;
	virtual void Serialize(Archive& ar) = 0;

public:
	//! recalculates the shape of the stiffness matrix
	bool CreateStiffness(bool breset);

public:
	FEM&		m_fem;	//!< reference the FE data structure

	// timers
	Timer	m_SolverTime;	//!< tracks time spent in solver

	// linear solver
	LinearSolver*		m_plinsolve;	//!< the linear solver

	// global stiffness matrix
	FEStiffnessMatrix*	m_pK;		//!< global stiffness matrix

	// convergence tolerances
	double	m_Rtol;			//!< residual tolerance
	double	m_Dtol;			//!< displacement tolerance
	double	m_Etol;			//!< energy tolerance
	double	m_Ptol;			//!< pressure tolerance
	double	m_LStol;		//!< line search tolerance
	double	m_LSmin;		//!< minimum line search step
	int		m_LSiter;		//!< max nr of line search iterations

	// BFGS parameters
	BFGSSolver	m_bfgs;		//!< BFGS solver parameters

	// counters
	int		m_nrhs;			//!< nr of right hand side evalutations
	int		m_niter;		//!< nr of quasi-newton iterations
	int		m_nref;			//!< nr of stiffness retormations
	int		m_naug;			//!< nr of augmentations
};

#endif // !defined(AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_)
