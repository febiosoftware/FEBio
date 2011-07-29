#pragma once

#include "FESolver.h"

//-----------------------------------------------------------------------------
//! The FELinearSolidSolver class solves linear (quasi-static) elasticity 
//! problems.
//!
class FELinearSolidSolver : public FESolver
{
public:
	//! constructor
	FELinearSolidSolver(FEM& fem);

	//! destructor (does nothing)
	~FELinearSolidSolver(){}

	//! Initialization
	bool Init();

	//! Solve the analysis step
	bool SolveStep(double time);

	//! Data serialization
	void Serialize(DumpFile& ar);

protected:
	//! calculate the residual
	void Residual();

	//! calculate the stiffness matrix
	bool StiffnessMatrix(); 

	//! form the stiffness matrix
	bool ReformStiffness();

	//! update solution
	void Update(vector<double>& u);

public:
	//! assemble element stiffness matrix
	void AssembleStiffness(matrix& ke, vector<int>& lm);

protected:
	vector<double>	m_u;	//!< nodal displacements
	vector<double>	m_R;	//!< right hand side
	vector<double>	m_d;	//!< prescribed displacements
};
