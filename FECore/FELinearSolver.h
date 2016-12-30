#pragma once
#include "FECore/FESolver.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEGlobalVector;
class FEGlobalMatrix;
class FELinearSystem;
class LinearSolver;

//-----------------------------------------------------------------------------
//! Abstract Base class for finite element solution algorithms that require the solution
//! of linear system of equations. This can be simple linear solution algorithm, but
//! also nonlinear algorithms that require solving a linear system of equations (e.g. Newton solvers)
class FELinearSolver : public FESolver
{
public:
	//! constructor
	FELinearSolver(FEModel* pfem);

	//! Set the degrees of freedom
	void SetDOF(vector<int>& dof);

	//! Get the number of equations
	int NumberOfEquations() const;

public: // from FESolver

	//! solve the step
	bool SolveStep(double time);

	//! Initialize and allocate data
	bool Init();

	//! Initialize equation numbers
	bool InitEquations();

	//! Clean up data
	void Clean();

	//! assemble global stiffness matrix
	void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

	//! Serialization
	void Serialize(DumpStream& ar);

public: // these functions need to be implemented by the derived class

	//! Evaluate the right-hand side "force" vector
	virtual void ForceVector(FEGlobalVector& R) = 0;

	//! Evaluate the stiffness matrix (new interface)
	virtual bool StiffnessMatrix(FELinearSystem& K) { return false; }

	//! Evaluate the stiffness matrix (old interface)
	virtual bool StiffnessMatrix();

	//! Update the model state
	virtual void Update(vector<double>& u);

protected: // some helper functions

	//! Reform the stiffness matrix
	bool ReformStiffness();

	//! Create and evaluate the stiffness matrix
	bool CreateStiffness();

protected:
	vector<double>		m_R;	//!< RHS vector
	vector<double>		m_u;	//!< vector containing prescribed values

private:
	LinearSolver*		m_pls;		//!< The linear equation solver
	FEGlobalMatrix*		m_pK;		//!< The global stiffness matrix
	int					m_neq;		//!< The number of equations (TODO: Get this from linear solver)

	vector<int>		m_dof;	//!< list of active degrees of freedom
	bool			m_breform;	//!< matrix reformation flag
};
