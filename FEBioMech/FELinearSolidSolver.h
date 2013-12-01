#pragma once
#include "FECore/FESolver.h"

class FEStiffnessMatrix;

//-----------------------------------------------------------------------------
//! The FELinearSolidSolver class solves linear (quasi-static) elasticity 
//! problems.
//!
class FELinearSolidSolver : public FESolver
{
public:
	//! constructor
	FELinearSolidSolver(FEModel* pfem);

	//! destructor
	~FELinearSolidSolver();

	//! Clean up
	virtual void Clean();

	//! Initialization
	bool Init();

	//! Solve the analysis step
	bool SolveStep(double time);

	//! Data serialization
	void Serialize(DumpFile& ar);

	//! Initialize linear equation system
	bool InitEquations();

protected:
	//! calculate the residual
	void Residual();

	//! recalculates the shape of the stiffness matrix
	bool CreateStiffness(bool breset);

	//! calculate the stiffness matrix
	bool StiffnessMatrix(); 

	//! form the stiffness matrix
	bool ReformStiffness();

	//! update solution
	void Update(vector<double>& u);

public:
	//! assemble element stiffness matrix
	void AssembleStiffness(matrix& ke, vector<int>& lm);

	//! assemble the element residual into the global residual
	void AssembleResidual(vector<int>& en, vector<int>& elm, vector<double>& fe, vector<double>& R);

private: // TODO: use this one instead
	void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke) { assert(false); }

public:
	double	m_Dtol;			//!< displacement tolerance

protected:
	vector<double>	m_u;	//!< nodal displacements
	vector<double>	m_R;	//!< right hand side
	vector<double>	m_d;	//!< prescribed displacements

public:
	LinearSolver*		m_plinsolve;	//!< the linear solver

	// global stiffness matrix
	FEStiffnessMatrix*	m_pK;		//!< global stiffness matrix
	int					m_neq;		//!< number of equations

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
