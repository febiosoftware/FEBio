// FESolver.h: interface for the FESolidSolver class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_)
#define AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/BFGSSolver.h"
#include "FECore/NonLinearSystem.h"
#include "FEStiffnessMatrix.h"
#include "fem.h"
#include "Timer.h"
#include "FECore/FEException.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! The solver class.

//! This class is responsible for solving a single time step in the FE
//! analysis, possibly using results from previous solutions.
//! In the future this class might become a base-class for different FE solvers.
//! For instances there could be a different solver for quasi-static, dynamic and
//! eigenvalue problems.
class FESolver : public NonLinearSystem
{
public:
	FESolver(FEM& fem);
	virtual ~FESolver();

	virtual bool Init();
	virtual bool SolveStep(double time) = 0;
	virtual void Serialize(DumpFile& ar) = 0;
	virtual void Clean();

public:
	// Initialize linear equation system (TODO: Is this the right place to do this?)
	virtual bool InitEquations() = 0;

	//! assemble the element residual into the global residual
	virtual void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R) {}

private:
	// These functions have to be overwritten from NonLinearSystem
	// but are not yet used.
	void Evaluate(vector<double>& R) { assert(false); }
	void Jacobian(SparseMatrix& K) { assert(false); }
	bool Converged() { assert(false); return false; }

public:
	//! recalculates the shape of the stiffness matrix
	bool CreateStiffness(bool breset);

	//! Get active domains
	int Domains() { return m_Dom.size(); }

	//! Get active domain
	FEDomain* Domain(int i) { return &(m_fem.m_mesh.Domain(m_Dom[i])); }

public:
	// FEM data
	FEM&		m_fem;	//!< reference the FE data structure

	// domain data
	vector<int>		m_Dom;	//!< list of active domains for this analysis

	// timers
	Timer	m_SolverTime;	//!< tracks time spent in solver

	// linear solver
	LinearSolver*		m_plinsolve;	//!< the linear solver

	// global stiffness matrix
	FEStiffnessMatrix*	m_pK;		//!< global stiffness matrix
	int					m_neq;		//!< number of equations

	// BFGS parameters
	BFGSSolver	m_bfgs;		//!< BFGS solver parameters

	// counters
	int		m_nrhs;			//!< nr of right hand side evalutations
	int		m_niter;		//!< nr of quasi-newton iterations
	int		m_nref;			//!< nr of stiffness retormations
	int		m_ntotref;
	int		m_naug;			//!< nr of augmentations
};

#endif // !defined(AFX_FESOLVER_H__EFCA10FB_7487_44A1_A81D_8B53BE3BECEA__INCLUDED_)
