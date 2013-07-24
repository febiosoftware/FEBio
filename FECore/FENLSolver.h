#pragma once
#include "BFGSSolver.h"
#include "FEParameterList.h"
#include "DumpFile.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
class FENLSolver : public FEParamContainer
{
public:
	//! constructor
	FENLSolver(FEModel& fem);

	//! destructor
	virtual ~FENLSolver();

	//! Get the FE model
	FEModel& GetFEModel();

public:
	virtual bool Init() = 0;

	virtual void Clean() = 0;

public:
	//! assemble into global residual (TODO: this is only used by rigid joints)
//	virtual void AssembleResidual(vector<int>& lm, vector<double>& fe, vector<double>& R) { assert(false); }

	//! assemble the element residual into the global residual
//	virtual void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R) = 0;

	//! assemble global stiffness matrix (TODO: this is only used by rigid joints)
	virtual void AssembleStiffness(vector<int>& elm, matrix& ke) { assert(false); }

	//! assemble global stiffness matrix
	virtual void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) = 0;

	// Initialize linear equation system (TODO: Is this the right place to do this?)
	virtual bool InitEquations() = 0;

	//! Solve an analysis step
	virtual bool SolveStep(double time) = 0;

public:
	//! Evaluate the state of the current system
	virtual void Evaluate(std::vector<double>& F) { assert(false); };

	//! Update the state of the sytem
	virtual void Update(std::vector<double>& u) { assert(false); };

protected:
	FEModel&	m_fem;

public: // TODO: temporary data that I would like to move elsewhere

	bool		m_bsymm;		//!< symmetry flag for linear solver allocation
	int			m_solvertype;	//!< defines the type of solver; 0=BFGs, 1-Hager-Zhang NLCG

	// counters
	int		m_nrhs;			//!< nr of right hand side evalutations
	int		m_niter;		//!< nr of quasi-newton iterations
	int		m_nref;			//!< nr of stiffness retormations
	int		m_ntotref;
	int		m_naug;			//!< nr of augmentations
};
