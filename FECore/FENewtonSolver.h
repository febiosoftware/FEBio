#pragma once
#include "FESolver.h"
#include "BFGSSolver.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;

//-----------------------------------------------------------------------------
//! This class defines the base class for Newton-type solvers. 
//! \todo This is a work in progress. I'd like to call this the BFGSSolver,
//!       but this name is already in use. This class should handle the common
//!       tasks in all solver classes that use the BFGSSolver.
class FENewtonSolver : public FESolver
{
public:
	//! constructor
	FENewtonSolver(FEModel* pfem);

	//! destrcutor
	~FENewtonSolver();

public: // overloaded from FESolver

	//! Initialization
	bool Init();

	//! Initialize linear equation system
	bool InitEquations();

	//! Clean up
	void Clean();

	//! serialization
	void Serialize(DumpFile& ar);

	//! Solve an analysis step
	bool SolveStep(double time);

protected:
	//! Performs a linesearch
	double LineSearch(double s);

	//! Do a Quasi-Newton step
	//! This is called from SolveStep and must be implemented by derived classes.
	virtual bool Quasin(double time) = 0;

public:	
	BFGSSolver	m_bfgs;			//!< BFGS solver parameters

public:
	// line search options
	double	m_LSmin;		//!< minimum line search step
	double	m_LStol;		//!< line search tolerance
	int		m_LSiter;		//!< max nr of line search iterations

	// linear solver data
	LinearSolver*		m_plinsolve;	//!< the linear solver
	int					m_neq;			//!< number of equations

	DECLARE_PARAMETER_LIST();
};
