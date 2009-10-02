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
	void Serialize(Archive& ar);

protected:
	//! calculate the residual
	void Residual();

	//! calculate the stiffness matrix
	bool StiffnessMatrix(); 

	//! form the stiffness matrix
	bool ReformStiffness();

	//! calculate the conductive element stiffness matrix
	void ConductionStiffness(FESolidElement& el, matrix& ke);

	//! calculate the capacitance element stiffness matrix
	void CapacitanceStiffness(FESolidElement& el, matrix& ke);

	//! update solution
	void Update();

	//! assemble element stiffness matrix
	void AssembleStiffness(matrix& ke, vector<int>& lm);

protected:
	vector<double>	m_T;	//!< temperature array
	vector<double>	m_Tp;	//!< previous temperatures
	vector<double>	m_R;	//!< residual
	vector<double>	m_u;	//!< prescribed temperatures
};
