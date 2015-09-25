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

protected:
	//! Performs a linesearch
	double LineSearch(double s);

public:	
	BFGSSolver	m_bfgs;			//!< BFGS solver parameters

public:
	// line search options
	double	m_LSmin;		//!< minimum line search step
	double	m_LStol;		//!< line search tolerance
	int		m_LSiter;		//!< max nr of line search iterations
};
