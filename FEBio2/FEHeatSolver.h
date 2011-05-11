#pragma once

#include "FESolver.h"

//-----------------------------------------------------------------------------
//! The FEHeatSolver solves heat-conduction problems
//! 

class FEHeatSolver : public FESolver
{
public:
	//! constructor
	FEHeatSolver(FEM& fem);

	//! one-time initialization
	bool Init();

	//! solve a timestep
	bool SolveStep(double time);

	//! serialize data
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

public:
	vector<double>	m_T;	//!< temperature array
	vector<double>	m_Tp;	//!< previous temperatures
	vector<double>	m_R;	//!< residual
	vector<double>	m_u;	//!< prescribed temperatures
};
