#pragma once
#include "NumCore/NonLinearSystem.h"
using namespace NumCore;

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
class FENLSolver : public NonLinearSystem
{
public:
	FENLSolver(FEModel& fem) : m_fem(fem) {}
	virtual ~FENLSolver() {}

public:
	//! assemble the element residual into the global residual
	virtual void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R) = 0;

	//! assemble an element stiffness into the residual
	virtual void AssembleResidual(vector<int>& lm, matrix& ke){};

	//! assemble global stiffness matrix
	virtual void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) = 0;

	FEModel& GetFEModel() { return m_fem; }

private:
	// TODO: I'm overriding these functions, but they are not used yet
	virtual void Evaluate(std::vector<double>& F) { assert(false); };
	virtual void Jacobian(SparseMatrix& K) { assert(false); };
	virtual void Update(std::vector<double>& u) { assert(false); };
	virtual bool Converged() { assert(false); return false; };

protected:
	FEModel&	m_fem;
};
