#pragma once
#include "NumCore/NonLinearSystem.h"
#include "NumCore/BFGSSolver.h"
#include "FEParameterList.h"
#include "DumpFile.h"
using namespace NumCore;

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
class FENLSolver : public NonLinearSystem, public FEParamContainer
{
public:
	FENLSolver(FEModel& fem) : m_fem(fem)
	{ 
		m_bsymm = true; // assume symmetric stiffness matrix
		m_solvertype = 0;
	}

	virtual ~FENLSolver() {}

	FEModel& GetFEModel() { return m_fem; }

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

private:
	// TODO: I'm overriding these functions, but they are not used yet
	virtual void Evaluate(std::vector<double>& F) { assert(false); };
	virtual void Jacobian(SparseMatrix& K) { assert(false); };
	virtual void Update(std::vector<double>& u) { assert(false); };
	virtual bool Converged() { assert(false); return false; };

protected:
	FEModel&	m_fem;

public: // TODO: temporary data that I would like to move elsewhere

	// BFGS parameters
	BFGSSolver	m_bfgs;		//!< BFGS solver parameters
	bool		m_bsymm;	//!< symmetry flag for linear solver allocation
	int			m_solvertype;	//!< defines the type of solver; 0=BFGs, 1-Hager-Zhang NLCG

	// counters
	int		m_nrhs;			//!< nr of right hand side evalutations
	int		m_niter;		//!< nr of quasi-newton iterations
	int		m_nref;			//!< nr of stiffness retormations
	int		m_ntotref;
	int		m_naug;			//!< nr of augmentations
};
