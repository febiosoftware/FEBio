#pragma once
#include "FECoreBase.h"
#include "DumpFile.h"
#include "Timer.h"
#include "matrix.h"
#include "vector.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! This is the base class for all FE solvers.

//! A class derived from FESolver implements a solver for a specific type
//! of physics problem. It takes the FEModel in its constructor and implements
//! the SolveStep function to solve the FE problem.
class FESolver : public FECoreBase
{
public:
	//! constructor
	FESolver(FEModel* pfem);

	//! destructor
	virtual ~FESolver();

	//! Get the FE model
	FEModel& GetFEModel();

public:
	//! Initialize solver data
	virtual bool Init();

	//! Can this be removed?
	virtual void Clean() = 0;

public:
	//! assemble global stiffness matrix (TODO: this is only used by rigid joints)
	virtual void AssembleStiffness(vector<int>& elm, matrix& ke) { assert(false); }

	//! assemble global stiffness matrix (TODO: this is only used by mortar contact)
	virtual void AssembleStiffness2(vector<int>& lmi, vector<int>& lmj, matrix& ke) { assert(false); }

	//! assemble global stiffness matrix
	virtual void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) = 0;

	// Initialize linear equation system (TODO: Is this the right place to do this?)
	// \todo Can I make this part of the Init function?
	virtual bool InitEquations() = 0;

	//! initialize the step
	virtual bool InitStep(double time);

	//! Solve an analysis step
	virtual bool SolveStep(double time) = 0;

public:
	//! Evaluate the state of the current system
	virtual void Evaluate(std::vector<double>& F) { assert(false); };

	//! Update the state of the sytem
	virtual void Update(std::vector<double>& u) { assert(false); };

protected:
	FEModel&	m_fem;

public: //TODO Move these parameters elsewhere
	bool		m_bsymm;		//!< symmetry flag for linear solver allocation

	// timers
	Timer	m_SolverTime;		//!< time spent in linear solver
	Timer	m_ReformTime;		//!< time reforming stiffness matrix
	Timer	m_StiffnessTime;	//!< time evaluating stiffness matrix
	Timer	m_RHSTime;			//!< time evaluating right-hand-side
	Timer	m_UpdateTime;		//!< time updating model

	// counters
	int		m_nrhs;			//!< nr of right hand side evalutations
	int		m_niter;		//!< nr of quasi-newton iterations
	int		m_nref;			//!< nr of stiffness retormations
	int		m_ntotref;		//!< nr of total stiffness reformations
	int		m_naug;			//!< nr of augmentations
};
