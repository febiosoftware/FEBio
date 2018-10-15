#pragma once
#include "FECoreBase.h"
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
class FECORE_API FESolver : public FECoreBase
{
	DECLARE_SUPER_CLASS(FESOLVER_ID);

public:
	//! constructor
	FESolver(FEModel* fem);

	//! destructor
	virtual ~FESolver();

public:
	//! Initialize solver data
	virtual bool Init() override;

	//! Data serialization
	void Serialize(DumpStream& ar) override;

	//! This is called by FEAnalaysis::Deactivate
	virtual void Clean();

	//! rewind the solver (This is called when the time step fails and needs to retry)
	virtual void Rewind() {}

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

	//! initialize the step (This is called before SolveStep)
	virtual bool InitStep(double time);

	//! Solve an analysis step
	virtual bool SolveStep() = 0;

	//! Update the state of the model
	virtual void Update(std::vector<double>& u) { assert(false); };

	//! Update the state of the model
	//! TODO: It might make more sense to move this to another class (e.g. FEAnalysis)
	virtual void UpdateModel() {}

	//! derived classes need to implement this.
	//! return true if the augmentations have converged
	virtual bool Augment() { return true; }

    //! Generate warnings if needed
    virtual void SolverWarnings() {}
    
public: //TODO Move these parameters elsewhere
	bool		m_bsymm;		//!< symmetry flag for linear solver allocation

	// counters
	int		m_nrhs;			//!< nr of right hand side evalutations
	int		m_niter;		//!< nr of quasi-newton iterations
	int		m_nref;			//!< nr of stiffness retormations
	int		m_ntotref;		//!< nr of total stiffness reformations

	// augmentation
	int		m_naug;			//!< nr of augmentations
	bool	m_baugment;		//!< do augmentations flag

	DECLARE_FECORE_CLASS();
};
